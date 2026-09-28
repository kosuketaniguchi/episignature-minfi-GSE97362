## =====================================================================
## 02_minfi_preprocess_qc.R
## 図2「DNAメチル化アレイの前処理とQC：IDATから解析用の行列まで」の minfi 実装
##
## 図2のブロックとの対応（[ ] 内が図の箱）
##   [IDAT + サンプルシート]  → read.metharray.sheet / read.metharray.exp
##   [検体QC]  検出率        → detectionP()（陰性コントロール法）
##             全体強度      → getQC()/plotQC()（M/U 中央値 log2）
##             変換率        → wateRmelon::bscon()
##             推定性別      → getSex()  ↔ サンプルシートの Sex
##             SNPプローブ   → getSnpBeta()（同一性・混入）
##             細胞組成      → FlowSorted.Blood.EPIC::estimateCellCounts2()（代替 EpiDISH）
##   [背景補正・色素補正]     → preprocessNoob(dyeCorr = TRUE)
##   [プローブのマスク]       → detectionP 不良、Zhou 2017 MASK_general（+集団別SNP）、
##                              dropLociWithSnps、交差反応、性染色体
##   [β/M 行列 + メタデータ + 解析記録] → getBeta / getM、pData、params.json、sessionInfo
##
## SeSAMe openSesame() との違い（README 参照）:
##   - 検出p は pOOBAH ではなく陰性コントロール法（minfi::detectionP）。pOOBAH 相当は ENmix::calcdetP(..., pt="oob")
##   - プローブマスクは自動適用されないので Zhou 2017 のマスク列を明示的に読む
##   - 色素補正は noob 内の dyeCorr（線形・参照式）。SeSAMe は非線形（dyeBiasNL）
## =====================================================================
suppressPackageStartupMessages({
  library(minfi)
  library(data.table)
  library(jsonlite)
})

## ---------------------------------------------------------------------
## 0. パラメータ（解析記録としてそのまま params.json に保存する）
## ---------------------------------------------------------------------
source("config.R")
params <- list(
  sample_sheet      = SAMPLESHEET,                             # 01 の出力
  out_dir           = PREPROC_DIR,
  array             = "450K",           # "450K" | "EPIC" | "EPICv2"
  genome            = "hg19",           # 450K/EPIC: hg19, EPICv2: hg38
  detP_method       = "m+u",            # minfi::detectionP の type
  detP_cutoff       = 0.01,
  sample_min_detected_frac = 0.95,      # 検体: 検出プローブ割合がこれ未満なら除外（=不良プローブ >5%）
  sample_badQC_cutoff      = 10.5,      # 検体: log2(M中央値), log2(U中央値) の平均がこれ未満なら除外（minfi 既定）
  sample_min_bscon         = 80,        # 検体: バイサルファイト変換率 (%)
  probe_max_failed_frac    = 0.05,      # プローブ: 検出p>cutoff の検体割合がこれ超なら除外
  noob_offset       = 15,
  noob_dyeMethod    = "single",         # "single"（1検体ごと）| "reference"
  apply_BMIQ        = FALSE,            # TRUE なら wateRmelon::BMIQ でプローブ型（I/II）補正
  mask_zhou_general = TRUE,             # Zhou et al. 2017 MASK_general
  mask_zhou_file    = file.path(ANNOT_DIR, "HM450.hg19.manifest.tsv.gz"),  # 00_download_data.R で取得
  mask_population_file = NULL,          # 例: file.path(ANNOT_DIR, "HM450.hg19.manifest.pop.tsv.gz")（EAS 列。日本人検体のとき）
  mask_population_col  = "MASK_snp5_EAS",
  mask_custom_probe_file = NULL,        # 日本人頻度（TogoVar 等）で作った自前リスト（1列, probeID）
  drop_snp_probes   = TRUE,             # minfi::dropLociWithSnps（CpG/SBE 位置の dbSNP、maf 0）
  drop_sex_chr      = TRUE,
  cell_composition  = TRUE,
  seed              = 20260928
)
## Zhou 2017 マニフェスト: https://zwdzwd.github.io/InfiniumAnnotation
##   450K: .../HM450/HM450.hg19.manifest.tsv.gz   EPIC: .../EPIC/EPIC.hg19.manifest.tsv.gz
##   EPICv2: .../EPICv2/EPICv2.hg38.manifest.tsv.gz  （集団別 SNP マスクは *.manifest.pop.tsv.gz、要確認）

dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(params$out_dir, "qc_plots"); dir.create(fig_dir, showWarnings = FALSE)
set.seed(params$seed)
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n")
## PDF 出力（日本語ラベル対応）：macOS は quartz、他は pdf(family = "Japan1")。pdf() の代わりにそのまま使える
if (capabilities("aqua")) quartzFonts(jp = quartzFont(rep("HiraginoSans-W3", 4)))
open_pdf <- function(file, width = 7, height = 7) {
  if (capabilities("aqua")) { quartz(file = file, type = "pdf", width = width, height = height, family = "jp"); par(family = "jp") }
  else pdf(file, width = width, height = height, family = "Japan1")
}

## ---------------------------------------------------------------------
## 1. IDAT + サンプルシート → RGChannelSet（図2 左端）
## ---------------------------------------------------------------------
targets <- fread(params$sample_sheet)
stopifnot("Basename" %in% colnames(targets))
targets[, Basename := file.path(IDAT_DIR, basename(Basename))]   # IDAT の場所は config.R から組み立てる
targets <- as.data.frame(targets)
log_line("IDAT 読込:", nrow(targets), "検体")
RGset <- read.metharray.exp(targets = targets, extended = FALSE, force = TRUE, verbose = FALSE)
sampleNames(RGset) <- targets$Sample_Name
## アレイ種別が自動認識されない場合（EPIC v2 など）は明示する
if (params$array == "EPICv2") annotation(RGset) <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
print(annotation(RGset))
ann <- getAnnotation(RGset)   # chr, pos, Name, Type(I/II), ...

## ---------------------------------------------------------------------
## 2. 検体QC（図2 中央上段）—— 自動処理できない部分。全部「記録」する
## ---------------------------------------------------------------------
qc <- data.table(Sample_Name = sampleNames(RGset))

## 2a. 検出率（陰性コントロール法）
detP <- detectionP(RGset, type = params$detP_method)
qc[, detected_frac := colMeans(detP < params$detP_cutoff)[Sample_Name]]
qc[, mean_detP     := colMeans(detP)[Sample_Name]]

## 2b. 全体強度（メチル化/非メチル化シグナル中央値）
MSet_raw <- preprocessRaw(RGset)
mu <- as.data.frame(getQC(MSet_raw))         # mMed, uMed (log2)
qc[, `:=`(mMed = mu$mMed[match(Sample_Name, rownames(mu))],
          uMed = mu$uMed[match(Sample_Name, rownames(mu))])]
qc[, intensity_ok := (mMed + uMed) / 2 >= params$sample_badQC_cutoff]
open_pdf(file.path(fig_dir, "01_plotQC_intensity.pdf")); plotQC(getQC(MSet_raw), badSampleCutoff = params$sample_badQC_cutoff); dev.off()

## 2c. バイサルファイト変換率（コントロールプローブ）
bs <- tryCatch(wateRmelon::bscon(RGset), error = function(e) { warning("bscon 失敗: ", conditionMessage(e)); rep(NA_real_, ncol(RGset)) })
qc[, bscon_pct := as.numeric(bs)[match(Sample_Name, sampleNames(RGset))]]
open_pdf(file.path(fig_dir, "02_control_probes.pdf"), width = 9, height = 6)
controlStripPlot(RGset, controls = c("BISULFITE CONVERSION I", "BISULFITE CONVERSION II"))
dev.off()

## 2d. 推定性別 ↔ 記録の照合（noob 後の X/Y 中央値強度から）
MSet_noob <- preprocessNoob(RGset, offset = params$noob_offset, dyeCorr = TRUE, dyeMethod = params$noob_dyeMethod)
GMset     <- mapToGenome(MSet_noob)
sexdf     <- as.data.frame(getSex(GMset, cutoff = -2))
qc[, `:=`(xMed = sexdf$xMed[match(Sample_Name, rownames(sexdf))],
          yMed = sexdf$yMed[match(Sample_Name, rownames(sexdf))],
          predictedSex = sexdf$predictedSex[match(Sample_Name, rownames(sexdf))])]
qc[, recordedSex := targets$Sex[match(Sample_Name, targets$Sample_Name)]]
qc[, sex_match := ifelse(is.na(recordedSex), NA, predictedSex == recordedSex)]
GMset <- addSex(GMset, getSex(GMset, cutoff = -2))          # plotSex() は addSex 済みの GenomicMethylSet を受け取る
open_pdf(file.path(fig_dir, "03_sex_check.pdf")); plotSex(GMset); dev.off()

## 2e. SNP プローブ（rs 番号プローブ）: 同一性（重複・検体入替）と混入
snpB <- getSnpBeta(RGset)                    # 450K: 65, EPIC v1: 59, EPIC v2: 65 程度
## 混入の目安: SNP β が 0.2–0.8 の「中間値」に落ちる割合。清浄な検体はヘテロ接合のみ中間になる
qc[, snp_intermediate_frac := colMeans(snpB > 0.2 & snpB < 0.8, na.rm = TRUE)[Sample_Name]]   # 参考値のみ（ヘテロ接合も数えるので混入の指標にはしない）
## 混入の指標：各 SNP の β が遺伝型の山（0 / 0.5 / 1）からどれだけ外れているかの平均。
## 清浄な検体は山に乗るので小さく、別人の DNA が混ざると 0.1〜0.35 や 0.65〜0.9 に値が散って大きくなる
snp_dev <- apply(snpB, 2, function(b) mean(pmin(abs(b), abs(b - 0.5), abs(b - 1)), na.rm = TRUE))
qc[, snp_genotype_dev := snp_dev[Sample_Name]]
## より厳密には ewastools::call_genotypes() → snp_outliers()（平均対数オッズ > -4 で混入疑い）
## 同一性: 検体間の SNP β 相関。同一人物・重複なら > 0.9 程度
snp_cor <- cor(snpB, use = "pairwise.complete.obs")
dup_pairs <- which(snp_cor > 0.9 & upper.tri(snp_cor), arr.ind = TRUE)
dup_dt <- data.table(sample1 = rownames(snp_cor)[dup_pairs[, 1]], sample2 = colnames(snp_cor)[dup_pairs[, 2]],
                     r = snp_cor[dup_pairs])
fwrite(dup_dt, file.path(params$out_dir, "qc_snp_duplicate_pairs.csv"))
open_pdf(file.path(fig_dir, "04_snp_probe_heatmap.pdf"), width = 8, height = 8)
pheatmap::pheatmap(snp_cor, show_rownames = ncol(snpB) < 80, show_colnames = FALSE, main = "SNP probe beta correlation")
dev.off()
## より厳密には ewastools::check_snp_agreement() / snp_outliers() を推奨

## 2f. 細胞組成（血液）: 共変量であり、極端な値は所見でもある
## 3 法を試し、成功したものをすべて保存する。qc には第一の成功法を採用し、複数あれば一致度を図示する。
##  A) FlowSorted.Blood.EPIC::estimateCellCounts2  IDOL 参照（Salas 2018/2022）。新しい版では referencePlatform は EPIC 固定で、
##     450K 入力は内部で共通プローブに揃えられる。450K には IDOLOptimizedCpGs450klegacy を渡す
##  B) minfi::estimateCellCounts + FlowSorted.Blood.450k  Houseman/Jaffe 法（450K ネイティブ参照。Gran = 顆粒球）
##  C) EpiDISH RPC + centDHSbloodDMC.m（333 DMC、7 分画。β 行列だけで動く）
if (params$cell_composition) {
  cc_list <- list()
  cc_list$IDOL <- tryCatch({
    suppressPackageStartupMessages(library(FlowSorted.Blood.EPIC))
    idol <- if (params$array == "450K") FlowSorted.Blood.EPIC::IDOLOptimizedCpGs450klegacy else FlowSorted.Blood.EPIC::IDOLOptimizedCpGs
    est <- estimateCellCounts2(RGset, compositeCellType = "Blood", processMethod = "preprocessNoob", probeSelect = "IDOL",
                               cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                               referencePlatform = "IlluminaHumanMethylationEPIC", IDOLOptimizedCpGs = idol)
    ## 版によって $prop（割合）または $counts で返る。空なら失敗扱い
    res <- if (!is.null(est$prop)) est$prop else if (!is.null(est$counts)) est$counts else est
    res <- as.data.frame(res)
    if (ncol(res) < 5 || nrow(res) != ncol(RGset)) stop("返り値の形が想定外: ", ncol(res), " 列 × ", nrow(res), " 行")
    res
  }, error = function(e) { warning("IDOL (estimateCellCounts2) 失敗: ", conditionMessage(e)); NULL })
  if (params$array == "450K") cc_list$Houseman450k <- tryCatch({
    if (!requireNamespace("FlowSorted.Blood.450k", quietly = TRUE)) BiocManager::install("FlowSorted.Blood.450k", update = FALSE, ask = FALSE)
    suppressPackageStartupMessages(library(FlowSorted.Blood.450k))
    est <- minfi::estimateCellCounts(RGset, compositeCellType = "Blood", referencePlatform = "IlluminaHumanMethylation450k",
                                     cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"))
    as.data.frame(est)
  }, error = function(e) { warning("Houseman/450k (estimateCellCounts) 失敗: ", conditionMessage(e)); NULL })
  cc_list$EpiDISH <- tryCatch({
    suppressPackageStartupMessages(library(EpiDISH))
    as.data.frame(epidish(beta.m = getBeta(MSet_noob), ref.m = EpiDISH::centDHSbloodDMC.m, method = "RPC")$estF)
  }, error = function(e) { warning("EpiDISH 失敗: ", conditionMessage(e)); NULL })
  cc_list <- Filter(Negate(is.null), cc_list)
  stopifnot(length(cc_list) > 0)

  ## 全法の結果を保存（列名に手法を付ける）
  cc_all <- Reduce(function(a, b) merge(a, b, by = "Sample_Name", all = TRUE),
                   lapply(names(cc_list), function(m) { d <- cc_list[[m]]; colnames(d) <- paste(m, colnames(d), sep = "."); d$Sample_Name <- rownames(cc_list[[m]]); d }))
  fwrite(cc_all, file.path(params$out_dir, "qc_cell_composition_all_methods.csv"))

  ## qc に採用するのは最初に成功した手法（IDOL > Houseman450k > EpiDISH）
  cc_method <- names(cc_list)[1]; cc <- cc_list[[1]]
  ## 手法間で列名を揃える（Neu / Gran / Neutro → Neu、EpiDISH の Eosino は Neu に合算して顆粒球扱い）
  harmonize <- function(d) {
    d <- as.data.frame(d); nm <- colnames(d)
    if ("Neutro" %in% nm) { d$Neu <- d$Neutro + if ("Eosino" %in% nm) d$Eosino else 0; d$Neutro <- NULL; d$Eosino <- NULL }
    if ("Gran" %in% nm) { d$Neu <- d$Gran; d$Gran <- NULL }
    if ("B" %in% nm) { d$Bcell <- d$B; d$B <- NULL }
    d[, intersect(c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"), colnames(d))]
  }
  cc_h <- lapply(cc_list, harmonize)
  cc <- cc_h[[1]]; cc$Sample_Name <- rownames(cc)
  ## 再実行時に古い細胞列が残っていると列名が衝突するので先に消す
  stale <- intersect(colnames(qc), c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu", "Gran", "B", "Neutro", "Eosino", "lymph_frac", "cell_method"))
  if (length(stale)) qc[, (stale) := NULL]
  qc <- merge(qc, cc, by = "Sample_Name", all.x = TRUE, sort = FALSE)
  qc[, cell_method := cc_method]
  qc[, lymph_frac := CD8T + CD4T + NK + Bcell]
  log_line("細胞組成: 採用 =", cc_method, "| 試行成功 =", paste(names(cc_list), collapse = ", "))

  ## 手法間の一致（2 法以上あるとき）：分画ごとの散布図と相関
  if (length(cc_h) >= 2) {
    open_pdf(file.path(fig_dir, "06_cell_composition_method_agreement.pdf"), width = 10, height = 7)
    par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
    m1 <- names(cc_h)[1]; m2 <- names(cc_h)[2]
    for (ct in intersect(colnames(cc_h[[1]]), colnames(cc_h[[2]]))) {
      x <- cc_h[[1]][[ct]]; y <- cc_h[[2]][rownames(cc_h[[1]]), ct]
      plot(x, y, pch = 19, cex = 0.6, col = adjustcolor("grey30", 0.6), xlab = m1, ylab = m2,
           main = sprintf("%s  r = %.2f", ct, cor(x, y, use = "complete.obs"))); abline(0, 1, lty = 2, col = "grey60")
    }
    dev.off()
  }
  ## 群別・年齢別の分布（所見としての細胞組成）
  grp <- targets$group[match(qc$Sample_Name, targets$Sample_Name)]; age <- targets$Age[match(qc$Sample_Name, targets$Sample_Name)]
  open_pdf(file.path(fig_dir, "07_cell_composition_by_group_age.pdf"), width = 11, height = 7)
  par(mfrow = c(2, 3), mar = c(7, 4, 2.5, 1))
  for (ct in c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")) if (ct %in% colnames(qc))
    boxplot(qc[[ct]] ~ grp, las = 2, cex.axis = 0.7, main = paste(ct, "(", cc_method, ")"), ylab = "推定割合", outline = TRUE)
  par(mfrow = c(2, 3), mar = c(4, 4, 2.5, 1))
  for (ct in c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")) if (ct %in% colnames(qc))
    plot(age, qc[[ct]], pch = 19, cex = 0.6, col = adjustcolor(ifelse(grp == "Control", "grey40", "#B9541F"), 0.6),
         xlab = "Age (years)", ylab = ct, main = paste(ct, "vs age（橙 = 症例）"))
  dev.off()
}

## 2g. 判定
qc[, fail_detection := detected_frac < params$sample_min_detected_frac]
qc[, fail_intensity := !intensity_ok]
qc[, fail_bscon     := !is.na(bscon_pct) & bscon_pct < params$sample_min_bscon]
qc[, flag_sex       := !is.na(sex_match) & !sex_match]
## 混入疑い：遺伝型の山からのずれが、コホートの中央値 + 5 MAD を超える検体（コホート内の相対基準）
## 旧基準「0.2〜0.8 の割合 > 0.5」はヘテロ接合の多い検体を誤って拾う（GSE97362 で 29 例が該当、全員が連続分布の上端だった）
qc[, flag_contam    := snp_genotype_dev > median(snp_genotype_dev, na.rm = TRUE) + 5 * mad(snp_genotype_dev, na.rm = TRUE)]
log_line(sprintf("SNP 遺伝型からのずれ: 中央値 %.3f, 最大 %.3f, 混入疑いの閾値 %.3f",
                 median(qc$snp_genotype_dev), max(qc$snp_genotype_dev),
                 median(qc$snp_genotype_dev) + 5 * mad(qc$snp_genotype_dev)))
## 細胞組成の極端値は「所見」として記録（除外はしない）：対照の分布から外れる検体（各分画で対照の中央値 ± 3 MAD 外）
if ("Neu" %in% colnames(qc)) {
  ctrl_idx <- targets$group[match(qc$Sample_Name, targets$Sample_Name)] == "Control"
  out_any <- rep(FALSE, nrow(qc))
  for (ct in intersect(c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"), colnames(qc))) {
    med <- median(qc[[ct]][ctrl_idx], na.rm = TRUE); madv <- mad(qc[[ct]][ctrl_idx], na.rm = TRUE)
    out_any <- out_any | (abs(qc[[ct]] - med) > 3 * madv)
  }
  qc[, flag_cell_extreme := out_any]
}
qc[, sample_pass    := !(fail_detection | fail_intensity | fail_bscon)]
## 性別不一致・混入疑いは「品質」ではなく「身元」の問題：自動除外せず、記録して個別に判断する
fwrite(qc, file.path(params$out_dir, "qc_samples.csv"))
log_line("検体QC: 合格", sum(qc$sample_pass), "/", nrow(qc),
         "| 性別不一致", sum(qc$flag_sex, na.rm = TRUE), "| 混入疑い", sum(qc$flag_contam, na.rm = TRUE))

keep_samples <- qc[sample_pass == TRUE & !(flag_sex %in% TRUE) & !(flag_contam %in% TRUE), Sample_Name]
RGset_f    <- RGset[, keep_samples]
MSet_noob  <- MSet_noob[, keep_samples]
detP       <- detP[, keep_samples]

## ---------------------------------------------------------------------
## 3. 背景補正・色素補正（noob）→ 任意で BMIQ（図2 中央下段）
## ---------------------------------------------------------------------
## noob は 2d で既に実施済み（MSet_noob）。ここで任意の BMIQ
if (params$apply_BMIQ) {
  beta_norm <- tryCatch(wateRmelon::BMIQ(MSet_noob),   # MethylSet メソッド：プローブ型 I/II の分布を揃えた β 行列
                        error = function(e) { warning("BMIQ 失敗（", conditionMessage(e), "）。noob の β を使います。"); getBeta(MSet_noob) })
} else {
  beta_norm <- getBeta(MSet_noob)
}
GRset <- ratioConvert(mapToGenome(MSet_noob), what = "both", keepCN = TRUE)  # GenomicRatioSet（M, CN も保持）

## ---------------------------------------------------------------------
## 4. プローブのマスク（図2 中央下段）
## ---------------------------------------------------------------------
probe_log <- list(start = nrow(GRset))
keep <- rep(TRUE, nrow(GRset)); names(keep) <- featureNames(GRset)

## 4a. 検出不良プローブ（検出p > cutoff の検体が probe_max_failed_frac 超）
failed_frac <- rowMeans(detP[featureNames(GRset), ] > params$detP_cutoff)
keep <- keep & (failed_frac <= params$probe_max_failed_frac)
probe_log$after_detP <- sum(keep)

## 4b. Zhou 2017 MASK_general（マッピング不良・反復配列・プローブ内 common SNP 等の推奨マスク）
if (params$mask_zhou_general && file.exists(params$mask_zhou_file)) {
  zm <- fread(params$mask_zhou_file, select = c("probeID", "MASK_general"))
  keep <- keep & !(names(keep) %in% zm[MASK_general == TRUE, probeID])
  probe_log$after_zhou_general <- sum(keep)
} else if (params$mask_zhou_general) warning("Zhou マニフェストが見つかりません: ", params$mask_zhou_file)

## 4c. 集団別 SNP マスク（東アジア集団 EAS）または自前リスト（TogoVar 等から作成）
if (!is.null(params$mask_population_file) && file.exists(params$mask_population_file)) {
  pm <- fread(params$mask_population_file)
  stopifnot(params$mask_population_col %in% colnames(pm))
  keep <- keep & !(names(keep) %in% pm[get(params$mask_population_col) == TRUE, probeID])
  probe_log$after_population_mask <- sum(keep)
}
if (!is.null(params$mask_custom_probe_file) && file.exists(params$mask_custom_probe_file)) {
  cm <- fread(params$mask_custom_probe_file, header = FALSE)[[1]]
  keep <- keep & !(names(keep) %in% cm)
  probe_log$after_custom_mask <- sum(keep)
}

## 4d. minfi 内蔵の dbSNP 注釈（CpG 位置・単一塩基伸長 SBE 位置の SNP）
if (params$drop_snp_probes) {
  snpinfo <- getSnpInfo(GRset)
  has_snp <- !is.na(snpinfo$CpG_rs) | !is.na(snpinfo$SBE_rs)
  keep <- keep & !has_snp
  probe_log$after_dropLociWithSnps <- sum(keep)
}

## 4e. 交差反応プローブ（Zhou の MASK_general に大半が含まれる。別途入れるなら maxprobes）
if (requireNamespace("maxprobes", quietly = TRUE)) {
  xr <- maxprobes::xreactive_probes(array_type = ifelse(params$array == "450K", "450K", "EPIC"))
  keep <- keep & !(names(keep) %in% xr)
  probe_log$after_crossreactive <- sum(keep)
}

## 4f. 性染色体
if (params$drop_sex_chr) {
  keep <- keep & !(as.character(seqnames(GRset)) %in% c("chrX", "chrY"))
  probe_log$after_sexchr <- sum(keep)
}
## 4g. EPIC v2 の重複プローブ（同一 CpG の _BC11/_TC21 等）は接尾辞を除いて平均する（Peters 2024 の推奨に準拠）
collapse_epicv2_replicates <- function(mat) {
  base_id <- sub("_[A-Z]{2}\\d{2}$", "", rownames(mat))
  if (!anyDuplicated(base_id)) return(mat)
  s <- rowsum(mat, group = base_id, na.rm = TRUE)          # 群ごとの合計
  n <- rowsum((!is.na(mat)) * 1, group = base_id)          # 群ごとの非NA数
  s / n                                                    # 平均
}

GRset_f   <- GRset[keep, ]
beta_norm <- beta_norm[featureNames(GRset_f), colnames(GRset_f)]
log_line("プローブフィルタ:", paste(names(probe_log), unlist(probe_log), sep = "=", collapse = " → "))

## ---------------------------------------------------------------------
## 5. 出力：β / M 行列、メタデータ、解析記録（図2 右端）
## ---------------------------------------------------------------------
beta <- beta_norm
M    <- log2(beta / (1 - beta))              # BMIQ 後も一貫させるため β から計算（β=getBeta の場合 getM と一致）
if (params$array == "EPICv2") { beta <- collapse_epicv2_replicates(beta); M <- log2(beta / (1 - beta)) }
pdat <- as.data.frame(colData(GRset_f))
pdat <- merge(pdat, qc, by = "Sample_Name", sort = FALSE)
probe_ann <- as.data.frame(ann[rownames(beta), intersect(c("chr", "pos", "strand", "Name", "Type", "Relation_to_Island", "UCSC_RefGene_Name"), colnames(ann))])

saveRDS(RGset_f,  file.path(params$out_dir, "RGset_filtered.rds"))
saveRDS(GRset_f,  file.path(params$out_dir, "GRset_noob_filtered.rds"))
fwrite(as.data.table(beta, keep.rownames = "probeID"), file.path(params$out_dir, "beta_matrix.csv.gz"))
fwrite(as.data.table(M,    keep.rownames = "probeID"), file.path(params$out_dir, "M_matrix.csv.gz"))
fwrite(pdat,      file.path(params$out_dir, "metadata_samples.csv"))
fwrite(as.data.table(probe_ann, keep.rownames = "probeID"), file.path(params$out_dir, "metadata_probes.csv.gz"))

record <- list(
  date = as.character(Sys.Date()), params = params, n_samples_in = nrow(targets), n_samples_out = ncol(beta),
  probe_filter_log = probe_log, array_annotation = as.list(annotation(RGset)),
  packages = sapply(c("minfi", "wateRmelon", "FlowSorted.Blood.EPIC", "EpiDISH", "limma"),
                    function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA))
)
write_json(record, file.path(params$out_dir, "params.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(params$out_dir, "sessionInfo.txt"))

## β 分布（プローブ型別）— 正規化の確認図
open_pdf(file.path(fig_dir, "05_beta_density.pdf"), width = 8, height = 5)
densityPlot(beta, sampGroups = if ("group" %in% colnames(pdat)) pdat$group else NULL, main = "noob (+BMIQ) beta")
dev.off()
log_line("完了:", nrow(beta), "CpG ×", ncol(beta), "検体 →", params$out_dir)
