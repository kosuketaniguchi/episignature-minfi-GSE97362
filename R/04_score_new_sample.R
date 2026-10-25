## =====================================================================
## 04_score_new_sample.R
## 「手元の 1 検体は Kabuki（KMT2D）か CHARGE（CHD7）か」を、03 で保存したシグネチャーで照合する。
##
## 流れ（図2 と同じ処理を新検体だけに適用 → 図3 の参照分布の上に置く）
##   IDAT 読込 → 最小限の検体QC（検出率・強度・推定性別・SNP 中間値・細胞組成）
##   → noob（検体内で完結する正規化なので参照コホートの再処理は不要）
##   → 2 本のシグネチャーでスコア（使えた CpG 数を報告。EPIC v1/v2 では一部欠ける）
##   → SVM 確率（任意）→ 参照分布（GSE97362 検証側）の上に描く
##
## 使い方：下の `new_basenames` に IDAT の Basename（_Grn/_Red.idat を除いたパス）を並べる。
##   デモ：GSE97362 の検体を「新患」に見立てる（demo_samples に Sample_Name を指定）
## =====================================================================
suppressPackageStartupMessages({ library(minfi); library(data.table); library(ggplot2) })
## 細胞組成は 02 と同じ IDOL を使う（FlowSorted.Blood.EPIC）。失敗時のみ EpiDISH

source("config.R")
sig_dir <- SIG_DIR
out_dir <- NEWSAMPLE_DIR; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 入力 -----------------------------------------------------------
demo_samples  <- c("KDM6A-1", "KMT2D-20")        # デモ：GSE97362 の検体名。実検体を使うときは NULL にする
new_basenames <- NULL                              # 例: c("/path/to/203xxxxx_R01C01")  ※ _Grn.idat を除く
array_type    <- "450K"                            # 新検体のアレイ: "450K" | "EPIC" | "EPICv2"
recorded_sex  <- NULL                              # 例: c("F")  記録があれば照合する

if (!is.null(demo_samples)) {
  ss <- fread(SAMPLESHEET)
  new_basenames <- file.path(IDAT_DIR, basename(ss[match(demo_samples, Sample_Name), Basename]))
  recorded_sex  <- ss[match(demo_samples, Sample_Name), Sex]
  new_names     <- paste0("DEMO_", demo_samples)
} else new_names <- basename(new_basenames)
stopifnot(length(new_basenames) > 0)

## ---- 1. 読込と最小限の QC ----------------------------------------------
RG <- read.metharray(new_basenames, force = TRUE)
sampleNames(RG) <- new_names
if (array_type == "EPICv2") annotation(RG) <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
detP <- detectionP(RG, type = "m+u")
qc   <- data.table(Sample_Name = new_names,
                   detected_frac = colMeans(detP < 0.01),
                   intensity     = rowMeans(as.matrix(as.data.frame(getQC(preprocessRaw(RG))))))
MS   <- preprocessNoob(RG, offset = 15, dyeCorr = TRUE, dyeMethod = "single")   # 02 と同じ設定
GM   <- mapToGenome(MS)
## getSex() は k-means を使うため 1 検体では動かない → 同じ計算（X/Y の log2 総強度中央値の差 < -2 で F）を自前で行う
sex_single <- function(GM, cutoff = -2) {
  CN  <- getCN(GM); ann <- getAnnotation(GM)
  xMed <- apply(CN[ann$chr == "chrX", , drop = FALSE], 2, median, na.rm = TRUE)
  yMed <- apply(CN[ann$chr == "chrY", , drop = FALSE], 2, median, na.rm = TRUE)
  data.frame(xMed = xMed, yMed = yMed, predictedSex = ifelse(yMed - xMed < cutoff, "F", "M"))
}
sx <- tryCatch(as.data.frame(getSex(GM, cutoff = -2)), error = function(e) sex_single(GM))
qc[, predictedSex := sx$predictedSex]
if (!is.null(recorded_sex)) qc[, sex_match := predictedSex == recorded_sex]
snpB <- getSnpBeta(RG)
## 混入の指標（02 と同じ）：SNP の β の、遺伝型の山（0 / 0.5 / 1）からの平均距離。02 の qc_samples.csv の分布と比べて判断する
qc[, snp_genotype_dev := apply(snpB, 2, function(b) mean(pmin(abs(b), abs(b - 0.5), abs(b - 1)), na.rm = TRUE))]
beta <- getBeta(MS)
if (array_type == "EPICv2") {                       # 重複プローブ（_BC11 等）を平均して 450K/EPIC v1 の ID に揃える
  base_id <- sub("_[A-Z]{2}\\d{2}$", "", rownames(beta))
  beta <- rowsum(beta, base_id, na.rm = TRUE) / rowsum((!is.na(beta)) * 1, base_id)
}
## 細胞組成：02 と同じ IDOL（estimateCellCounts2）を優先。参照データと一緒に処理する方式なので 1 検体でも動く。
## 失敗時のみ EpiDISH（02 では比較用）。列名は 02 と同じ CD8T/CD4T/NK/Bcell/Mono/Neu に揃える
cc <- tryCatch({
  suppressPackageStartupMessages(library(FlowSorted.Blood.EPIC))
  idol <- if (array_type == "450K") FlowSorted.Blood.EPIC::IDOLOptimizedCpGs450klegacy else FlowSorted.Blood.EPIC::IDOLOptimizedCpGs
  est <- estimateCellCounts2(RG, compositeCellType = "Blood", processMethod = "preprocessNoob", probeSelect = "IDOL",
                             cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                             referencePlatform = "IlluminaHumanMethylationEPIC", IDOLOptimizedCpGs = idol)
  r <- as.data.frame(if (!is.null(est$prop)) est$prop else est$counts)
  if (ncol(r) < 5) stop("返り値が空"); attr(r, "method") <- "IDOL"; r
}, error = function(e) {
  warning("IDOL 失敗（", conditionMessage(e), "）。EpiDISH で代替します。02 の参照値（IDOL）とは推定法が異なる点に注意。")
  r <- as.data.frame(EpiDISH::epidish(beta.m = beta, ref.m = EpiDISH::centDHSbloodDMC.m, method = "RPC")$estF)
  r <- data.frame(CD8T = r$CD8T, CD4T = r$CD4T, NK = r$NK, Bcell = r$B, Mono = r$Mono, Neu = r$Neutro + r$Eosino, row.names = rownames(r))
  attr(r, "method") <- "EpiDISH"; r
})
qc <- cbind(qc, as.data.table(round(cc[qc$Sample_Name, c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")], 3)))
qc[, `:=`(lymph_frac = CD8T + CD4T + NK + Bcell, cell_method = attr(cc, "method"))]
qc[, qc_pass := detected_frac >= 0.95 & intensity >= 10.5]
print(qc)

## ---- 2. シグネチャー照合 ------------------------------------------------
sigs <- list(KMT2D = readRDS(file.path(sig_dir, "signature_KMT2D.rds")), CHD7 = readRDS(file.path(sig_dir, "signature_CHD7.rds")))
score_one <- function(b, sig) {
  ids <- intersect(sig$cpgs$probeID, rownames(b)); s <- sig$cpgs[match(ids, probeID)]
  z <- (b[ids, , drop = FALSE] - s$ctrl_mean) / pmax(s$ctrl_sd, 0.01)
  sc <- colMeans(z * sign(s$delta_beta), na.rm = TRUE)
  ## SVM は全 CpG が必要：欠けた CpG は発見対照の平均で埋める（中立値）ことを明示して補う
  p_svm <- if (!is.null(sig$svm)) {
    X <- matrix(sig$cpgs$ctrl_mean, nrow = ncol(b), ncol = nrow(sig$cpgs), byrow = TRUE, dimnames = list(colnames(b), sig$cpgs$probeID))
    X[, ids] <- t(b[ids, , drop = FALSE])
    attr(predict(sig$svm, X, probability = TRUE), "probabilities")[, "case"]
  } else NA_real_
  data.table(Sample_Name = colnames(b), signature = sig$name, score = sc,
             cutoff = sig$threshold$cutoff, positive = sc > sig$threshold$cutoff,
             z_vs_heldout_ctrl = (sc - sig$threshold$heldout_ctrl_mean) / sig$threshold$heldout_ctrl_sd,
             svm_p_case = p_svm, cpg_used = length(ids), cpg_total = nrow(sig$cpgs))
}
res <- rbindlist(lapply(sigs, function(s) score_one(beta, s)))
print(res[, .(Sample_Name, signature, score = round(score, 2), cutoff = round(cutoff, 2), positive,
              z = round(z_vs_heldout_ctrl, 1), svm_p = round(svm_p_case, 3), cpg = paste0(cpg_used, "/", cpg_total))])
fwrite(res, file.path(out_dir, "new_sample_scores.csv")); fwrite(qc, file.path(out_dir, "new_sample_qc.csv"))

## ---- 3. 参照分布の上に置く -------------------------------------------------
ref <- fread(file.path(sig_dir, "reference_scores.csv"))[used_for == "none"]   # 選抜に使っていない 165 例
grp_order <- c("Control", "CHARGE_CHD7_LOF", "CHARGE_CHD7_P_LP", "CHARGE_CHD7_VUS", "KS1_KMT2D_VUS", "KS1_KMT2D_P_LP", "KS1_KMT2D_LOF", "KS2_KDM6A")
pal <- c(Control = "#BDBDBD", CHARGE_CHD7_LOF = "#2E5E7E", CHARGE_CHD7_P_LP = "#4F86A8", CHARGE_CHD7_VUS = "#9DBBD0",
         KS1_KMT2D_VUS = "#E3B08F", KS1_KMT2D_P_LP = "#D0783C", KS1_KMT2D_LOF = "#B9541F", KS2_KDM6A = "#7B2D8E")
ref[, group := factor(group, levels = grp_order)]
## デモ検体は参照側からは除いて描く（同じ点を二度描かない）
if (!is.null(demo_samples)) ref <- ref[!Sample_Name %in% demo_samples]
newpts <- dcast(res, Sample_Name ~ signature, value.var = "score")
p <- ggplot() +
  geom_vline(xintercept = sigs$KMT2D$threshold$cutoff, linetype = 2, colour = "grey40") +
  geom_hline(yintercept = sigs$CHD7$threshold$cutoff,  linetype = 2, colour = "grey40") +
  geom_point(data = ref, aes(score_KMT2D, score_CHD7, colour = group), size = 1.9, alpha = 0.7) +
  scale_colour_manual(values = pal, drop = FALSE, name = "参照（GSE97362 検証側）") +
  geom_point(data = newpts, aes(KMT2D, CHD7), shape = 23, size = 4.5, fill = "#FFD400", colour = "black", stroke = 1) +
  geom_text(data = newpts, aes(KMT2D, CHD7, label = Sample_Name), nudge_y = 0.12, size = 3, fontface = "bold") +
  labs(x = "KMT2D（Kabuki）シグネチャースコア", y = "CHD7（CHARGE）シグネチャースコア",
       title = "新検体の照合：参照分布の上に置く", subtitle = "◆ = 新検体。破線 = 未使用対照の平均 + 3SD") +
  theme_bw(base_size = 11)
if (capabilities("aqua")) {
  quartzFonts(jp = quartzFont(rep("HiraginoSans-W3", 4)))
  quartz(file = file.path(out_dir, "new_sample_on_reference.pdf"), type = "pdf", width = 9, height = 5.5, family = "jp")
  print(p + theme(text = element_text(family = "jp")))
} else {
  pdf(file.path(out_dir, "new_sample_on_reference.pdf"), width = 9, height = 5.5, family = "Japan1")
  print(p + theme(text = element_text(family = "Japan1")))
}
dev.off()
cat("完了 →", out_dir, "\n")
