## =====================================================================
## 03_build_episignatures.R
## GSE97362（02 の出力）から KMT2D（Kabuki）と CHD7（CHARGE）の 2 本のエピシグネチャーを構築し、
## 新しい検体の照合に使える形（04 が読む）で保存する。あわせて原稿の図3を描く。
##
## 入力: PREPROC_DIR/{beta_matrix.csv.gz, M_matrix.csv.gz, metadata_samples.csv}（config.R）
## 出力: SIG_DIR/（config.R）
##   signature_KMT2D.rds / signature_CHD7.rds    照合に必要な一式（下記 build_signature の返り値）
##   reference_scores.csv                       全 235 検体の 2 本のスコア（04 で参照分布として描く）
##   fig3_two_signature_scores.pdf              検証セットの 2 次元スコア図（横 KMT2D、縦 CHD7）
##   fig3b_score_by_group_{KMT2D,CHD7}.pdf      群別スコア（LOO 併記）
##   check_celltype_overlap.csv / check_score_vs_lymph.pdf   細胞組成との切り分け
##
## 発見セット（GEO の "sample type" で定義済み）:
##   KMT2D: KS1_KMT2D_LOF 11 例 + Control 11 例（年齢を揃えた対照）
##   CHD7 : CHARGE_CHD7_LOF 19 例 + Control 29 例
## 検証セット = それぞれの発見セットに使っていない検体（相手側の発見コホートも含む）
## =====================================================================
suppressPackageStartupMessages({ library(data.table); library(limma); library(ggplot2); library(pheatmap) })
set.seed(20260928)

source("config.R")
in_dir  <- PREPROC_DIR
out_dir <- SIG_DIR; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

beta <- as.matrix(fread(file.path(in_dir, "beta_matrix.csv.gz")), rownames = "probeID")
M    <- as.matrix(fread(file.path(in_dir, "M_matrix.csv.gz")),    rownames = "probeID")
pd   <- fread(file.path(in_dir, "metadata_samples.csv"))
stopifnot(all(pd$Sample_Name == colnames(beta)))
pd[is.na(Age), Age := median(pd[group == "Control", Age], na.rm = TRUE)]
cell_cols <- intersect(c("CD8T", "CD4T", "NK", "Bcell", "Mono"), colnames(pd))   # Neu は参照（合計 1）として外す

## ---------------------------------------------------------------------
## 選抜パラメータ（params としてシグネチャーに同梱する）
##  主解析: 共変量は年齢・性別のみ。発見コホートは年齢を揃えてあり、n=22 に細胞 5 分画を足すと自由度を失う上、
##          疾患で本当に変わる分画（Kabuki の Bcell）を共変量にすると疾患効果を吸い取る（媒介変数の調整）。
##          細胞組成の影響は事後確認（細胞型 CpG との重なり、スコア vs リンパ球割合）で担保する。
##  感度分析: 細胞 5 分画も共変量に入れた版を別に作り、主解析との CpG の重なりを報告する。
## ---------------------------------------------------------------------
sel <- list(adj_p = 0.05, min_delta_beta = 0.05, r_max = 0.85, n_max = 200,
            covariates = c("Age", "Sex"), stat = "limma on M-values", score = "mean signed z vs discovery controls")
sel_cells <- modifyList(sel, list(covariates = c("Age", "Sex", cell_cols)))
## 図出力：日本語ラベルのため、macOS は quartz、それ以外は pdf(family="Japan1") を使う
## （既定の pdf() は日本語が "..." に化け、cairo_pdf は XQuartz が無いと "failed to load cairo DLL" になる）
if (capabilities("aqua")) quartzFonts(jp = quartzFont(rep("HiraginoSans-W3", 4)))   # 古い macOS は "HiraKakuProN-W3"
jp_family <- if (capabilities("aqua")) "jp" else "Japan1"
jp_theme  <- theme(text = element_text(family = jp_family))
open_pdf <- function(file, width, height) {
  if (capabilities("aqua")) { quartz(file = file, type = "pdf", width = width, height = height, family = "jp"); par(family = "jp") }
  else pdf(file, width = width, height = height, family = "Japan1")
}
gg_save <- function(file, plot, width, height) { open_pdf(file, width, height); print(plot + jp_theme); dev.off() }

prune_correlated <- function(ids, mat, r_max, n_max) {
  keep <- character(0)
  for (p in ids) {
    if (length(keep) >= n_max) break
    if (!length(keep)) { keep <- p; next }
    r <- suppressWarnings(cor(mat[p, ], t(mat[keep, , drop = FALSE])))
    if (all(abs(r) <= r_max, na.rm = TRUE)) keep <- c(keep, p)
  }
  keep
}
score_samples <- function(b, sig) {   # b: β行列（行 = CpG）。sig$cpgs にある CpG のうち b に存在するものだけ使う
  ids <- intersect(sig$cpgs$probeID, rownames(b)); s <- sig$cpgs[match(ids, probeID)]
  z <- (b[ids, , drop = FALSE] - s$ctrl_mean) / pmax(s$ctrl_sd, 0.01)
  list(score = colMeans(z * sign(s$delta_beta), na.rm = TRUE), n_used = length(ids), n_total = nrow(sig$cpgs))
}

## ---------------------------------------------------------------------
## シグネチャー構築（発見セットのみを使う）
## ---------------------------------------------------------------------
build_signature <- function(name, cases, ctrls, sel) {
  disc <- c(cases, ctrls)
  pdd  <- pd[match(disc, Sample_Name)]
  pdd[, case := factor(fifelse(Sample_Name %in% cases, "case", "ctrl"), levels = c("ctrl", "case"))]
  cov  <- sel$covariates[sapply(sel$covariates, function(v) length(unique(na.omit(pdd[[v]]))) > 1)]
  design <- model.matrix(as.formula(paste("~ case +", paste(cov, collapse = " + "))), data = pdd)
  stopifnot(nrow(design) == length(disc))
  fit <- eBayes(lmFit(M[, disc], design))
  tt  <- as.data.table(topTable(fit, coef = "casecase", number = Inf, sort.by = "none"), keep.rownames = "probeID")
  tt[, `:=`(delta_beta = rowMeans(beta[, cases, drop = FALSE]) - rowMeans(beta[, ctrls, drop = FALSE]),
            ctrl_mean  = rowMeans(beta[, ctrls, drop = FALSE]),
            ctrl_sd    = apply(beta[, ctrls, drop = FALSE], 1, sd))]
  ## 閾値の感度：候補数のグリッドを表示（どこで選抜が崩れるかを見る）
  grid <- CJ(adj_p = c(0.01, 0.05), min_db = c(0.05, 0.10))[, n := mapply(function(a, d) tt[adj.P.Val < a & abs(delta_beta) >= d, .N], adj_p, min_db)]
  cat(sprintf("[%s] 候補数グリッド: %s\n", name, paste(sprintf("p<%.2f&|Δβ|≥%.2f: %d", grid$adj_p, grid$min_db, grid$n), collapse = " | ")))
  cand <- tt[adj.P.Val < sel$adj_p & abs(delta_beta) >= sel$min_delta_beta][order(-abs(t))]
  ids  <- prune_correlated(cand$probeID, beta[, disc, drop = FALSE], sel$r_max, sel$n_max)
  cpgs <- cand[probeID %in% ids]
  if (nrow(cpgs) < 5) warning(sprintf("[%s] シグネチャー CpG が %d 個しかありません。閾値か n を見直してください。", name, nrow(cpgs)))
  ## SVM（線形）：発見セットで学習。04 で新検体に適用する（CpG が 5 個未満なら作らない）
  svm_fit <- if (requireNamespace("e1071", quietly = TRUE) && nrow(cpgs) >= 5)
    e1071::svm(t(beta[cpgs$probeID, disc, drop = FALSE]), pdd$case, kernel = "linear", probability = TRUE, cost = 1) else NULL
  sig <- list(name = name, gene = name, cases = cases, ctrls = ctrls, cpgs = cpgs, svm = svm_fit, params = sel,
              n_candidates = nrow(cand), built = as.character(Sys.Date()), platform = "IlluminaHumanMethylation450k (hg19)")
  ## 選抜に使っていない対照でスコア分布 → 陽性閾値（平均 + 3SD）
  held_ctrl <- setdiff(pd[group == "Control", Sample_Name], disc)
  sc <- score_samples(beta[, held_ctrl, drop = FALSE], sig)$score
  sig$threshold <- list(heldout_ctrl_mean = mean(sc), heldout_ctrl_sd = sd(sc), cutoff = mean(sc) + 3 * sd(sc), n_heldout_ctrl = length(sc))
  cat(sprintf("[%s] 候補 %d → シグネチャー %d CpG（低メチル化 %d / 高メチル化 %d）; 陽性閾値 %.2f（未使用対照 %d 例）\n",
              name, nrow(cand), nrow(cpgs), sum(cpgs$delta_beta < 0), sum(cpgs$delta_beta > 0), sig$threshold$cutoff, length(sc)))
  sig
}

## 共変量の判断（2026-09-28、対照内の「スコア vs リンパ球割合」の相関で決めた）:
##   KMT2D: 年齢・性別のみ。細胞共変量なしで r = 0.08（交絡なし）。入れると n=22 と Bcell の共線性で 200 → 8 CpG に崩れる
##   CHD7 : 細胞 5 分画も入れる。入れないと対照内で r = −0.69（好中球の多い検体でスコアが上がる）、入れると 0.15。感度 16/16 は不変
sig_KMT2D <- build_signature("KMT2D",
  cases = pd[set == "discovery_KMT2D" & group == "KS1_KMT2D_LOF", Sample_Name],
  ctrls = pd[set == "discovery_KMT2D" & group == "Control", Sample_Name], sel)
sig_CHD7 <- build_signature("CHD7",
  cases = pd[set == "discovery_CHD7" & group == "CHARGE_CHD7_LOF", Sample_Name],
  ctrls = pd[set == "discovery_CHD7" & group == "Control", Sample_Name], sel_cells)
saveRDS(sig_KMT2D, file.path(out_dir, "signature_KMT2D.rds"))
saveRDS(sig_CHD7,  file.path(out_dir, "signature_CHD7.rds"))
fwrite(sig_KMT2D$cpgs, file.path(out_dir, "signature_KMT2D_cpgs.csv"))
fwrite(sig_CHD7$cpgs,  file.path(out_dir, "signature_CHD7_cpgs.csv"))
cat("2 本のシグネチャーで重なる CpG:", length(intersect(sig_KMT2D$cpgs$probeID, sig_CHD7$cpgs$probeID)), "\n")

## 感度分析：共変量の扱いを逆にした版（保存はするが 04 では使わない）
cat("\n--- 感度分析（共変量を逆にした版：KMT2D は細胞あり、CHD7 は細胞なし）---\n")
sig_KMT2D_alt <- build_signature("KMT2D", sig_KMT2D$cases, sig_KMT2D$ctrls, sel_cells)
sig_CHD7_alt  <- build_signature("CHD7",  sig_CHD7$cases,  sig_CHD7$ctrls,  sel)
sens <- data.table(signature = c("KMT2D", "CHD7"),
                   covariates_primary = c("Age+Sex", "Age+Sex+cells"),
                   n_primary = c(nrow(sig_KMT2D$cpgs), nrow(sig_CHD7$cpgs)),
                   n_alt = c(nrow(sig_KMT2D_alt$cpgs), nrow(sig_CHD7_alt$cpgs)),
                   overlap = c(length(intersect(sig_KMT2D$cpgs$probeID, sig_KMT2D_alt$cpgs$probeID)),
                               length(intersect(sig_CHD7$cpgs$probeID,  sig_CHD7_alt$cpgs$probeID))))
print(sens); fwrite(sens, file.path(out_dir, "sensitivity_cell_covariates.csv"))
saveRDS(list(KMT2D = sig_KMT2D_alt, CHD7 = sig_CHD7_alt), file.path(out_dir, "signatures_alternative_covariates.rds"))

## ---------------------------------------------------------------------
## 全検体のスコア（参照分布として保存。発見セットの検体は "used_for" に印）
## ---------------------------------------------------------------------
ref <- pd[, .(Sample_Name, group, set, Age, Sex)]
ref[, score_KMT2D := score_samples(beta, sig_KMT2D)$score[Sample_Name]]
ref[, score_CHD7  := score_samples(beta, sig_CHD7)$score[Sample_Name]]
ref[, used_for := fcase(Sample_Name %in% c(sig_KMT2D$cases, sig_KMT2D$ctrls), "KMT2D_discovery",
                        Sample_Name %in% c(sig_CHD7$cases,  sig_CHD7$ctrls),  "CHD7_discovery", default = "none")]
ref[, pos_KMT2D := score_KMT2D > sig_KMT2D$threshold$cutoff]
ref[, pos_CHD7  := score_CHD7  > sig_CHD7$threshold$cutoff]
if ("lymph_frac" %in% colnames(pd)) ref[, lymph_frac := pd$lymph_frac[match(Sample_Name, pd$Sample_Name)]]

## LOO（発見セットの検体は 1 例抜いて再選抜したスコアで評価する）
loo_score <- function(sig, s) {
  cases <- setdiff(sig$cases, s); ctrls <- setdiff(sig$ctrls, s); disc <- c(cases, ctrls)
  pdd <- pd[match(disc, Sample_Name)]; pdd[, case := factor(fifelse(Sample_Name %in% cases, "case", "ctrl"), levels = c("ctrl", "case"))]
  cov <- sig$params$covariates[sapply(sig$params$covariates, function(v) length(unique(na.omit(pdd[[v]]))) > 1)]
  des <- model.matrix(as.formula(paste("~ case +", paste(cov, collapse = " + "))), data = pdd)
  tt  <- as.data.table(topTable(eBayes(lmFit(M[, disc], des)), coef = "casecase", number = Inf, sort.by = "none"), keep.rownames = "probeID")
  tt[, `:=`(delta_beta = rowMeans(beta[, cases, drop = FALSE]) - rowMeans(beta[, ctrls, drop = FALSE]),
            ctrl_mean = rowMeans(beta[, ctrls, drop = FALSE]), ctrl_sd = apply(beta[, ctrls, drop = FALSE], 1, sd))]
  cand <- tt[adj.P.Val < sig$params$adj_p & abs(delta_beta) >= sig$params$min_delta_beta][order(-abs(t))]
  ids  <- prune_correlated(cand$probeID, beta[, disc, drop = FALSE], sig$params$r_max, sig$params$n_max)
  if (!length(ids)) return(NA_real_)
  score_samples(beta[, s, drop = FALSE], list(cpgs = cand[probeID %in% ids]))$score
}
ref[used_for == "KMT2D_discovery", score_KMT2D_loo := sapply(Sample_Name, function(s) loo_score(sig_KMT2D, s))]
ref[used_for == "CHD7_discovery",  score_CHD7_loo  := sapply(Sample_Name, function(s) loo_score(sig_CHD7,  s))]
fwrite(ref, file.path(out_dir, "reference_scores.csv"))

## 陽性率（検証側のみ）
pos <- ref[used_for == "none", .(n = .N, pos_KMT2D = sum(pos_KMT2D), pos_CHD7 = sum(pos_CHD7)), by = group][order(group)]
print(pos); fwrite(pos, file.path(out_dir, "positives_by_group.csv"))

## ---------------------------------------------------------------------
## 図3：2 次元スコア図（検証側の検体。発見セットは LOO スコアで薄く重ねる）
## ---------------------------------------------------------------------
grp_order <- c("Control", "CHARGE_CHD7_LOF", "CHARGE_CHD7_P_LP", "CHARGE_CHD7_VUS",
               "KS1_KMT2D_VUS", "KS1_KMT2D_P_LP", "KS1_KMT2D_LOF", "KS2_KDM6A")
pal <- c(Control = "#9E9E9E", CHARGE_CHD7_LOF = "#2E5E7E", CHARGE_CHD7_P_LP = "#4F86A8", CHARGE_CHD7_VUS = "#9DBBD0",
         KS1_KMT2D_VUS = "#E3B08F", KS1_KMT2D_P_LP = "#D0783C", KS1_KMT2D_LOF = "#B9541F", KS2_KDM6A = "#7B2D8E")
plot_dt <- copy(ref)
plot_dt[used_for == "KMT2D_discovery", score_KMT2D := score_KMT2D_loo]   # 発見側は LOO 値で表示（循環を避ける）
plot_dt[used_for == "CHD7_discovery",  score_CHD7  := score_CHD7_loo]
plot_dt[, group := factor(group, levels = grp_order)]
plot_dt[, shown_as := fifelse(used_for == "none", "検証（選抜に不使用）", "発見セット（LOO スコア）")]
p <- ggplot(plot_dt, aes(score_KMT2D, score_CHD7, colour = group, shape = shown_as)) +
  geom_vline(xintercept = sig_KMT2D$threshold$cutoff, linetype = 2, colour = "grey40") +
  geom_hline(yintercept = sig_CHD7$threshold$cutoff,  linetype = 2, colour = "grey40") +
  geom_point(size = 2.3, alpha = 0.85) +
  scale_colour_manual(values = pal, drop = FALSE) + scale_shape_manual(values = c(16, 1)) +
  labs(x = "KMT2D（Kabuki）シグネチャースコア", y = "CHD7（CHARGE）シグネチャースコア", colour = NULL, shape = NULL,
       title = "GSE97362：2 本のエピシグネチャーによる照合（minfi, noob）",
       subtitle = "破線：選抜に使っていない対照の平均 + 3SD") +
  theme_bw(base_size = 11) + theme(legend.position = "right")
gg_save(file.path(out_dir, "fig3_two_signature_scores.pdf"), p, 9, 5.5)

for (g in c("KMT2D", "CHD7")) {
  sc <- paste0("score_", g); dd <- plot_dt[!is.na(get(sc))]
  q <- ggplot(dd, aes(group, .data[[sc]], colour = group, shape = shown_as)) +
    geom_hline(yintercept = get(paste0("sig_", g))$threshold$cutoff, linetype = 2, colour = "grey40") +
    geom_jitter(width = 0.15, size = 2, alpha = 0.85) + scale_colour_manual(values = pal, drop = FALSE) + scale_shape_manual(values = c(16, 1)) +
    labs(x = NULL, y = paste(g, "シグネチャースコア"), colour = NULL, shape = NULL) +
    theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
  gg_save(file.path(out_dir, sprintf("fig3b_score_by_group_%s.pdf", g)), q, 7, 4)
}

## ---------------------------------------------------------------------
## 事後確認：シグネチャーは細胞組成の差を拾っていないか
##  (1) 細胞型識別 CpG（EpiDISH の 333 DMC、IDOL ライブラリ）との重なり
##  (2) スコア vs リンパ球割合（検証側の対照で相関が強ければ危険信号）
## ---------------------------------------------------------------------
ct_sets <- list()
if (requireNamespace("EpiDISH", quietly = TRUE)) ct_sets$EpiDISH_333DMC <- rownames(EpiDISH::centDHSbloodDMC.m)
if (requireNamespace("FlowSorted.Blood.EPIC", quietly = TRUE)) {
  ct_sets$IDOL_450klegacy <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs450klegacy
  ct_sets$IDOL_EPIC <- FlowSorted.Blood.EPIC::IDOLOptimizedCpGs
}
ov <- rbindlist(lapply(names(ct_sets), function(n) data.table(
  celltype_set = n, n_set = length(ct_sets[[n]]),
  overlap_KMT2D = length(intersect(sig_KMT2D$cpgs$probeID, ct_sets[[n]])),
  overlap_CHD7  = length(intersect(sig_CHD7$cpgs$probeID,  ct_sets[[n]])))))
print(ov); fwrite(ov, file.path(out_dir, "check_celltype_overlap.csv"))
## 主解析と逆版で、対照内の「スコア vs リンパ球割合」相関と陽性率を比較する（共変量の判断根拠）
if ("lymph_frac" %in% colnames(ref)) {
  d <- ref[used_for == "none"]
  sigs_cmp <- list(KMT2D_primary = sig_KMT2D, KMT2D_alt = sig_KMT2D_alt, CHD7_primary = sig_CHD7, CHD7_alt = sig_CHD7_alt)
  cmp <- rbindlist(lapply(names(sigs_cmp), function(n) {
    s <- sigs_cmp[[n]]; sc <- score_samples(beta[, d$Sample_Name], s)$score; ctrl <- d$group == "Control"
    gene <- sub("_.*", "", n); tgt <- grepl(gene, d$group) & grepl("LOF|P_LP", d$group)
    data.table(signature = n, covariates = paste(s$params$covariates, collapse = "+"), n_cpg = nrow(s$cpgs),
               r_ctrl_lymph = round(cor(sc[ctrl], d$lymph_frac[ctrl]), 2),
               pos_ctrl = sum(sc[ctrl] > s$threshold$cutoff), n_ctrl = sum(ctrl),
               pos_LOF_PLP = sum(sc[tgt] > s$threshold$cutoff), n_LOF_PLP = sum(tgt))
  }))
  print(cmp); fwrite(cmp, file.path(out_dir, "check_cell_confounding_by_version.csv"))
  open_pdf(file.path(out_dir, "check_score_vs_lymph.pdf"), 9, 4.2); par(mfrow = c(1, 2))
  for (g in c("KMT2D", "CHD7")) {
    d <- ref[used_for == "none"]; sc <- d[[paste0("score_", g)]]
    r_ctrl <- cor(sc[d$group == "Control"], d$lymph_frac[d$group == "Control"], use = "complete.obs")
    plot(d$lymph_frac, sc, pch = 19, cex = 0.7, col = adjustcolor(pal[as.character(d$group)], 0.8),
         xlab = "推定リンパ球割合 (CD8T+CD4T+NK+Bcell)", ylab = paste(g, "score"), main = sprintf("%s: 対照内の r = %.2f", g, r_ctrl))
    abline(h = get(paste0("sig_", g))$threshold$cutoff, lty = 2, col = "grey40")
  }
  dev.off()
}
cat("完了 →", out_dir, "\n")
