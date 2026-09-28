## =====================================================================
## 05_heatmap_episignatures.R
## 03 で作った KMT2D / CHD7 シグネチャーの CpG をヒートマップで示す。
##
## A. 主図（heatmap_validation.pdf）：CpG の選抜に使っていない検体だけを表示する（QC 通過例のみ）
##   行 = シグネチャー CpG（KMT2D 低/高メチル化、CHD7 低/高メチル化の 4 ブロック。ブロック内はクラスタリング）
##   列 = 群ごとに区切り、群内はその群に対応するシグネチャーのスコア順
##   値 = 各 CpG の「発見コホート対照の平均からのずれ」z = (β − 対照平均) / 対照SD（±3 で打ち切り）
##        青 = 対照より低い、赤 = 高い、灰 = 対照並み。value = "beta" で β そのもの（単色の濃淡）に切替
##   上段 = 2 本のシグネチャースコア（破線 = 陽性閾値。閾値超えは黒、未満は灰）
## 補足図（heatmap_discovery_circular.pdf）：選抜に使った検体で描いた図。必ず綺麗に分かれるが検証にはならない
## B. 教師なし階層クラスタリング（heatmap_cluster_{KMT2D,CHD7,both}.pdf）：シグネチャーごとに 1 枚。
##   列（検体）・行（CpG）とも群ラベルを使わずにクラスタリングし、群は色帯として後から重ねるだけ。
##   樹形図を 2 つに切ったクラスタと群の対応表（cluster_vs_group_*.csv）も出す
## C. 図3A 用の小さなヒートマップ（heatmap_panelA_KMT2D.pdf）：KMT2D シグネチャーのみ、対照は無作為 12 例、
##   KMT2D P/LP・VUS・KDM6A は選抜に使っていない全例。幅 約 90 mm。抽出した検体は heatmap_panelA_samples.csv
##   C だけ作り直すときは、冒頭（入力の読込と rows の作成）と B の save_both_devices() の定義を実行してから C を実行する
##
## 入力: PREPROC_DIR/beta_matrix.csv.gz、SIG_DIR/{signature_*.rds, reference_scores.csv}（config.R）
## 出力: SIG_DIR/heatmap_*.pdf / .png、heatmap_*_matrix.csv.gz、cluster_vs_group_*.csv
## 文字は英語（日本語フォントの問題を避けるため。図の説明は本文・凡例で行う）
## =====================================================================
suppressPackageStartupMessages(library(data.table))
for (p in c("ComplexHeatmap", "circlize"))
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(grid) })
ht_opt$message <- FALSE

source("config.R")
sig_dir <- SIG_DIR
value   <- "z"      # "z"（対照からのずれ、推奨）| "beta"
z_clip  <- 3

beta <- as.matrix(fread(file.path(PREPROC_DIR, "beta_matrix.csv.gz")), rownames = "probeID")
ref  <- fread(file.path(sig_dir, "reference_scores.csv"))
sigs <- list(KMT2D = readRDS(file.path(sig_dir, "signature_KMT2D.rds")),
             CHD7  = readRDS(file.path(sig_dir, "signature_CHD7.rds")))
thr  <- sapply(sigs, function(s) s$threshold$cutoff)

## ---- 行：シグネチャー CpG（2 本のシグネチャーで同じ CpG があっても別の行として扱う）----
rows <- rbindlist(lapply(names(sigs), function(g) {
  s <- sigs[[g]]$cpgs
  data.table(row_id = paste(g, s$probeID, sep = ":"), probeID = s$probeID, signature = g,
             direction = fifelse(s$delta_beta < 0, "hypo", "hyper"),
             ctrl_mean = s$ctrl_mean, ctrl_sd = s$ctrl_sd)
}))
rows <- rows[probeID %in% rownames(beta)]
rows[, block := factor(paste(signature, direction),
                       levels = c("KMT2D hypo", "KMT2D hyper", "CHD7 hypo", "CHD7 hyper"))]
cat("CpG 数:", paste(names(table(rows$block)), table(rows$block), sep = " = ", collapse = ", "), "\n")

make_matrix <- function(samples) {
  b <- beta[rows$probeID, samples, drop = FALSE]
  m <- if (value == "z") pmin(pmax((b - rows$ctrl_mean) / pmax(rows$ctrl_sd, 0.01), -z_clip), z_clip) else b
  rownames(m) <- rows$row_id; m
}
## 色：z は発散型（青 − 灰 − 赤、中央は無彩色）、β は単色の濃淡
col_fun <- if (value == "z") colorRamp2(c(-z_clip, 0, z_clip), c("#2a78d6", "#f0efec", "#e34948")) else
                             colorRamp2(c(0, 1), c("#f0efec", "#1F4E8C"))
legend_title <- if (value == "z") "z vs discovery controls" else "beta"

## 群の色（家系 = 色相、P/LP と VUS = 濃淡。配色チェッカーで確認済み）と表示名
grp_col <- c(Control = "#9E9E9E",
             CHARGE_CHD7_LOF = "#1F5A8C", CHARGE_CHD7_P_LP = "#1F5A8C", CHARGE_CHD7_VUS = "#3F8FDA",
             KS1_KMT2D_LOF = "#A8431A", KS1_KMT2D_P_LP = "#A8431A", KS1_KMT2D_VUS = "#E07F45",
             KS2_KDM6A = "#7B2D8E")
grp_lab <- c(Control = "Control", CHARGE_CHD7_LOF = "CHD7 LOF", CHARGE_CHD7_P_LP = "CHD7 P/LP", CHARGE_CHD7_VUS = "CHD7 VUS",
             KS1_KMT2D_LOF = "KMT2D LOF", KS1_KMT2D_P_LP = "KMT2D P/LP", KS1_KMT2D_VUS = "KMT2D VUS", KS2_KDM6A = "KDM6A")

## ---- ヒートマップ本体 ----
draw_heatmap <- function(cols, split, title, file_stub, show_scores = TRUE, width = 10, height = 7) {
  split <- droplevels(split)                       # 空の群があっても落ちないように
  m <- make_matrix(cols$Sample_Name)
  fwrite(as.data.table(m, keep.rownames = "row_id"), file.path(sig_dir, paste0(file_stub, "_matrix.csv.gz")))
  top <- if (show_scores) {
    pt <- function(g) {
      x <- cols[[paste0("score_", g)]]
      anno_points(x, ylim = range(c(x, thr[g]), na.rm = TRUE) + c(-0.3, 0.3), size = unit(1.1, "mm"), height = unit(13, "mm"),
                  gp = gpar(col = ifelse(x > thr[g], "black", "grey65")), axis_param = list(gp = gpar(fontsize = 6)))
    }
    HeatmapAnnotation(Group = as.character(cols$group), `KMT2D score` = pt("KMT2D"), `CHD7 score` = pt("CHD7"),
                      col = list(Group = grp_col), show_legend = FALSE, gap = unit(1.5, "mm"),
                      annotation_name_gp = gpar(fontsize = 7), annotation_name_side = "left", annotation_name_rot = 0)
  } else HeatmapAnnotation(Group = as.character(cols$group), col = list(Group = grp_col), show_legend = FALSE,
                           annotation_name_gp = gpar(fontsize = 7), annotation_name_side = "left", annotation_name_rot = 0)
  ht <- Heatmap(m, name = legend_title, col = col_fun, top_annotation = top,
                row_split = rows$block, cluster_row_slices = FALSE, cluster_rows = TRUE, show_row_dend = FALSE,
                row_title_gp = gpar(fontsize = 8), row_title_rot = 0,
                column_split = split, cluster_columns = FALSE, cluster_column_slices = FALSE,
                column_title_gp = gpar(fontsize = 8), column_title_rot = 90,
                show_row_names = FALSE, show_column_names = FALSE,
                column_gap = unit(1.2, "mm"), row_gap = unit(1.2, "mm"), border = FALSE,
                heatmap_legend_param = list(direction = "horizontal", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7),
                                            legend_width = unit(35, "mm")))
  lgd_pts <- if (show_scores) list(Legend(labels = c("above threshold (control mean + 3SD)", "below threshold"), type = "points", pch = 16,
                                          legend_gp = gpar(col = c("black", "grey65")), title = "Score", title_gp = gpar(fontsize = 8),
                                          labels_gp = gpar(fontsize = 7), nrow = 1)) else list()
  render <- function() {
    draw(ht, column_title = title, column_title_gp = gpar(fontsize = 10, fontface = "bold"),
         heatmap_legend_side = "bottom", annotation_legend_list = lgd_pts, annotation_legend_side = "bottom", merge_legend = TRUE)
    if (show_scores) for (g in c("KMT2D", "CHD7")) for (i in seq_len(nlevels(split)))
      decorate_annotation(paste(g, "score"), slice = i,
        grid.lines(c(0, 1), unit(c(thr[g], thr[g]), "native"), gp = gpar(lty = 2, col = "grey30", lwd = 0.8)))
  }
  pdf(file.path(sig_dir, paste0(file_stub, ".pdf")), width = width, height = height); render(); dev.off()
  png(file.path(sig_dir, paste0(file_stub, ".png")), width = width, height = height, units = "in", res = 200); render(); dev.off()
  cat("保存:", file.path(sig_dir, paste0(file_stub, ".pdf")), "\n")
}

## ---- 主図：選抜に使っていない検体 ----
val_levels <- c("Control", "CHARGE_CHD7_P_LP", "CHARGE_CHD7_VUS", "KS1_KMT2D_P_LP", "KS1_KMT2D_VUS", "KS2_KDM6A")
vd <- ref[used_for == "none" & group %in% val_levels]
vd[, group := factor(group, levels = val_levels)]
## 群内の並び：CHARGE 群は CHD7 スコア、Kabuki 群は KMT2D スコア、対照は大きい方のスコアの降順
vd[, ord := fcase(grepl("CHD7", group), -score_CHD7,
                  grepl("KMT2D|KDM6A", group), -score_KMT2D,
                  default = -pmax(score_KMT2D, score_CHD7))]
setorder(vd, group, ord)
split_v <- factor(grp_lab[as.character(vd$group)], levels = grp_lab[val_levels])
draw_heatmap(vd, split_v, sprintf("Episignature CpGs in %d samples not used for CpG selection (GSE97362)", nrow(vd)),
             "heatmap_validation")

## ---- 補足図：選抜に使った検体（循環。検証にはならないことを示すため）----
dd <- rbind(ref[used_for == "KMT2D_discovery"][, panel := fifelse(group == "Control", "KMT2D disc.: control", "KMT2D disc.: LOF")],
            ref[used_for == "CHD7_discovery"][,  panel := fifelse(group == "Control", "CHD7 disc.: control",  "CHD7 disc.: LOF")])
dd[, panel := factor(panel, levels = c("KMT2D disc.: control", "KMT2D disc.: LOF", "CHD7 disc.: control", "CHD7 disc.: LOF"))]
setorder(dd, panel, -score_KMT2D)
draw_heatmap(dd, dd$panel, "Discovery samples used for CpG selection (circular: separation here is not validation)",
             "heatmap_discovery_circular", show_scores = FALSE, width = 7, height = 7)


## =====================================================================
## B. 教師なし階層クラスタリング：シグネチャーごとに 1 枚
##   距離 = ユークリッド距離、連結法 = Ward.D2（EpiSign 系の論文でよく使われる組み合わせ）
##   検体 = そのシグネチャーの CpG 選抜に使っていない検体すべて
##     KMT2D の図には CHD7 の発見コホート（CHARGE LOF と対照）も入る → 他疾患を拾わないかの確認
##     CHD7 の図には KMT2D の発見コホート（Kabuki LOF と対照）も入る
##   樹形図を k = 2 で切り、スコア平均の高い方を "Signature-like"、低い方を "Control-like" と呼ぶ
## =====================================================================
clust_dist <- "euclidean"; clust_method <- "ward.D2"; k_clusters <- 2
disp_group <- function(g) fcase(g == "Control", "Control",
  g %in% c("CHARGE_CHD7_LOF", "CHARGE_CHD7_P_LP"), "CHD7 LOF/P/LP", g == "CHARGE_CHD7_VUS", "CHD7 VUS",
  g %in% c("KS1_KMT2D_LOF", "KS1_KMT2D_P_LP"), "KMT2D LOF/P/LP", g == "KS1_KMT2D_VUS", "KMT2D VUS",
  g == "KS2_KDM6A", "KDM6A", default = "Other")
disp_levels <- c("Control", "CHD7 LOF/P/LP", "CHD7 VUS", "KMT2D LOF/P/LP", "KMT2D VUS", "KDM6A")
disp_col  <- c(Control = "#9E9E9E", `CHD7 LOF/P/LP` = "#1F5A8C", `CHD7 VUS` = "#3F8FDA",
               `KMT2D LOF/P/LP` = "#A8431A", `KMT2D VUS` = "#E07F45", KDM6A = "#7B2D8E")
clust_col <- c(`Signature-like` = "#333333", `Control-like` = "#D9D9D9")

cluster_heatmap <- function(g) {
  smp <- ref[used_for != paste0(g, "_discovery")]
  rr  <- rows[signature == g]
  b   <- beta[rr$probeID, smp$Sample_Name, drop = FALSE]
  m   <- if (value == "z") pmin(pmax((b - rr$ctrl_mean) / pmax(rr$ctrl_sd, 0.01), -z_clip), z_clip) else b
  rownames(m) <- rr$probeID
  hc  <- hclust(dist(t(m), method = clust_dist), method = clust_method)     # 群ラベルは使わない
  cl  <- cutree(hc, k = k_clusters)
  sc  <- smp[[paste0("score_", g)]]
  sig_cl <- as.integer(names(which.max(tapply(sc, cl, mean))))
  cl_lab <- ifelse(cl == sig_cl, "Signature-like", "Control-like")
  grp <- factor(disp_group(smp$group), levels = disp_levels)
  tab <- as.data.frame.matrix(table(grp, factor(cl_lab, levels = names(clust_col))))
  cat(sprintf("\n[%s] 教師なしクラスタ（k = %d）× 群\n", g, k_clusters)); print(tab)
  fwrite(as.data.table(tab, keep.rownames = "group"), file.path(sig_dir, sprintf("cluster_vs_group_%s.csv", g)))
  ## 注釈名は図ごとに一意にし（閾値線を描くため）、表示は annotation_label で揃える
  annos <- setNames(list(grp, cl_lab,
    anno_points(sc, ylim = range(c(sc, thr[g]), na.rm = TRUE) + c(-0.3, 0.3), size = unit(1, "mm"), height = unit(13, "mm"),
                gp = gpar(col = ifelse(sc > thr[g], "black", "grey65")), axis_param = list(gp = gpar(fontsize = 6)))),
    paste0(g, c("_group", "_cluster", "_score")))
  top <- do.call(HeatmapAnnotation, c(annos, list(
    col = setNames(list(disp_col, clust_col), paste0(g, c("_group", "_cluster"))),
    annotation_label = c("Group", "Unsupervised cluster", paste(g, "score")),
    annotation_name_gp = gpar(fontsize = 7), annotation_name_side = "left", annotation_name_rot = 0,
    annotation_legend_param = setNames(list(
      list(title = "Group", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7), nrow = 1),
      list(title = "Unsupervised cluster", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7), nrow = 1)),
      paste0(g, c("_group", "_cluster"))),
    gap = unit(1.2, "mm"))))
  Heatmap(m, name = paste0(g, "_values"), col = col_fun, top_annotation = top,
          cluster_columns = hc, column_dend_height = unit(18, "mm"),
          cluster_rows = TRUE, clustering_distance_rows = clust_dist, clustering_method_rows = clust_method,
          show_row_dend = TRUE, row_dend_width = unit(8, "mm"),
          show_row_names = FALSE, show_column_names = FALSE,
          column_title = sprintf("%s signature (%d CpGs) x %d samples not used for its CpG selection (unsupervised: %s, %s)",
                                 g, nrow(m), ncol(m), clust_dist, clust_method),
          column_title_gp = gpar(fontsize = 9, fontface = "bold"),
          heatmap_legend_param = list(title = legend_title, direction = "horizontal", title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 7), legend_width = unit(35, "mm")))
}
hts <- list(KMT2D = cluster_heatmap("KMT2D"), CHD7 = cluster_heatmap("CHD7"))

lgd_score <- list(Legend(labels = c("above threshold (control mean + 3SD)", "below threshold"), type = "points", pch = 16,
                         legend_gp = gpar(col = c("black", "grey65")), title = "Score", title_gp = gpar(fontsize = 8),
                         labels_gp = gpar(fontsize = 7), nrow = 1))
draw_one <- function(g, newpage = TRUE) {
  draw(hts[[g]], newpage = newpage, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
       annotation_legend_list = lgd_score, merge_legend = TRUE)
  decorate_annotation(paste0(g, "_score"),
    grid.lines(c(0, 1), unit(c(thr[g], thr[g]), "native"), gp = gpar(lty = 2, col = "grey30", lwd = 0.8)))
}
save_both_devices <- function(stub, width, height, render) {
  pdf(file.path(sig_dir, paste0(stub, ".pdf")), width = width, height = height); render(); dev.off()
  png(file.path(sig_dir, paste0(stub, ".png")), width = width, height = height, units = "in", res = 200); render(); dev.off()
  cat("保存:", file.path(sig_dir, paste0(stub, ".pdf")), "\n")
}
for (g in names(hts)) save_both_devices(paste0("heatmap_cluster_", g), 10, 6.5, function() draw_one(g))
## 2 枚を上下に並べた 1 ページ（A: KMT2D、B: CHD7）
save_both_devices("heatmap_cluster_both", 10, 12.5, function() {
  grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 1)))
  for (i in 1:2) {
    pushViewport(viewport(layout.pos.row = i)); draw_one(names(hts)[i], newpage = FALSE)
    grid.text(LETTERS[i], x = unit(3, "mm"), y = unit(1, "npc") - unit(3, "mm"), just = c("left", "top"),
              gp = gpar(fontsize = 14, fontface = "bold"))
    upViewport()
  }
})

## =====================================================================
## C. 図3A 用の小さなヒートマップ（KMT2D シグネチャーのみ、列を絞る）
##   「シグネチャーとはこういう模様」を見せるパネル。B パネル（04 の図：参照分布の上の◆）と並べて使う。
##   列：対照（選抜に使っていない対照から無作為に n_ctrl 例）、KMT2D P/LP・KMT2D VUS・KDM6A（選抜に使っていない全例）
##       群内は KMT2D スコアの降順。上段にスコアと閾値（破線）、B パネルで◆にした検体に◆印
##   行：KMT2D シグネチャー CpG を「KMT2D LOF で低い／高い」の 2 ブロックに分け、ブロック内はクラスタリング
##   大きさ：幅 約 90 mm（誌面の片段幅）で文字 6〜7 pt になるように作る
## =====================================================================
panelA <- list(
  n_ctrl   = 12,                                                   # 対照の数（無作為）
  seed     = 20260928,                                             # 無作為抽出の乱数シード（凡例に書く）
  groups   = c("KS1_KMT2D_P_LP", "KS1_KMT2D_VUS", "KS2_KDM6A"),   # "CHARGE_CHD7_P_LP" を足すと特異性も示せる
  mark     = c("KDM6A-1", "KMT2D-20"),                             # 04 のデモ検体（B パネルの◆）。不要なら character(0)
  width    = 3.6, height = 4.2,                                    # インチ（3.6 in ≒ 91 mm）
  footnote = TRUE)                                                 # 図の下に対照の抽出方法を小さく書く

set.seed(panelA$seed)
ctrl_pool <- ref[used_for == "none" & group == "Control", Sample_Name]
ctrl_pick <- sort(sample(ctrl_pool, min(panelA$n_ctrl, length(ctrl_pool))))
pa <- ref[used_for == "none" & (Sample_Name %in% ctrl_pick | group %in% panelA$groups)]
pa_levels <- c("Control", panelA$groups)
pa[, group := factor(group, levels = pa_levels)]
setorder(pa, group, -score_KMT2D)
pa[, positive := score_KMT2D > thr["KMT2D"]]
fwrite(pa[, .(Sample_Name, group, score_KMT2D, positive, marked_in_B = Sample_Name %in% panelA$mark)],
       file.path(sig_dir, "heatmap_panelA_samples.csv"))
cat(sprintf("\n[図3A] 対照 %d / %d 例を無作為抽出（seed %d）。群別の陽性数:\n", length(ctrl_pick), length(ctrl_pool), panelA$seed))
print(pa[, .(n = .N, positive = sum(positive)), by = group])

## 値の行列（KMT2D シグネチャーのみ）
rrA <- rows[signature == "KMT2D"]
bA  <- beta[rrA$probeID, pa$Sample_Name, drop = FALSE]
mA  <- if (value == "z") pmin(pmax((bA - rrA$ctrl_mean) / pmax(rrA$ctrl_sd, 0.01), -z_clip), z_clip) else bA
rownames(mA) <- rrA$probeID
lab_lo <- sprintf("Lower in\nKMT2D LOF\n(%d CpGs)",  sum(rrA$direction == "hypo"))
lab_hi <- sprintf("Higher in\nKMT2D LOF\n(%d CpGs)", sum(rrA$direction == "hyper"))
row_blk   <- factor(ifelse(rrA$direction == "hypo", lab_lo, lab_hi), levels = c(lab_lo, lab_hi))
pa_lab    <- c(Control = "Control", KS1_KMT2D_P_LP = "KMT2D\nP/LP", KS1_KMT2D_VUS = "KMT2D\nVUS",
               KS2_KDM6A = "KDM6A", CHARGE_CHD7_P_LP = "CHD7\nP/LP")
col_split <- droplevels(factor(pa_lab[as.character(pa$group)], levels = pa_lab[pa_levels]))

## 上段の注釈：◆印（任意）→ 群の色帯 → KMT2D スコア
scA <- pa$score_KMT2D
is_mark <- pa$Sample_Name %in% panelA$mark
annA <- list(); labA <- character(0)
if (length(panelA$mark)) {
  annA$A_mark <- anno_simple(ifelse(is_mark, "mark", "none"), col = c(mark = "white", none = "white"), border = FALSE,
                             pch = ifelse(is_mark, 23, NA), pt_gp = gpar(fill = "#FFD400", col = "black"),
                             pt_size = unit(2.2, "mm"), height = unit(3.2, "mm"))
  labA <- c(labA, "Shown in B")
}
annA$A_group <- as.character(pa$group); labA <- c(labA, "Group")
annA$A_score <- anno_points(scA, ylim = range(c(scA, thr["KMT2D"]), na.rm = TRUE) + c(-0.3, 0.3),
                            size = unit(1.3, "mm"), height = unit(10, "mm"),
                            gp = gpar(col = ifelse(scA > thr["KMT2D"], "black", "grey65")),
                            axis_param = list(gp = gpar(fontsize = 6)))
labA <- c(labA, "KMT2D score")
topA <- do.call(HeatmapAnnotation, c(annA, list(
  col = list(A_group = grp_col), show_legend = FALSE, annotation_label = labA,
  simple_anno_size = unit(2.8, "mm"), gap = unit(1, "mm"),
  annotation_name_gp = gpar(fontsize = 7), annotation_name_side = "left", annotation_name_rot = 0)))

htA <- Heatmap(mA, name = "panelA", col = col_fun, top_annotation = topA,
               row_split = row_blk, cluster_row_slices = FALSE, cluster_rows = TRUE, show_row_dend = FALSE,
               row_title_gp = gpar(fontsize = 7), row_title_rot = 0,
               column_split = col_split, cluster_columns = FALSE, cluster_column_slices = FALSE,
               column_title_gp = gpar(fontsize = 7), column_title_rot = 0,
               show_row_names = FALSE, show_column_names = FALSE,
               column_gap = unit(1, "mm"), row_gap = unit(1, "mm"),
               heatmap_legend_param = list(title = legend_title, direction = "horizontal", title_position = "topcenter",
                                           title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6),
                                           legend_width = unit(28, "mm")))
lgdA <- list(Legend(labels = c("above threshold", "below"), type = "points", pch = 16, size = unit(1.8, "mm"),
                    legend_gp = gpar(col = c("black", "grey65")), title = "KMT2D score", title_position = "topcenter",
                    title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6), nrow = 1))
renderA <- function() {
  ## padding = 下・左・上・右。右は KDM6A の列見出しのはみ出し分、下は注記の分
  draw(htA, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", annotation_legend_list = lgdA,
       merge_legend = TRUE, padding = unit(c(if (panelA$footnote) 9 else 2, 2, 2, 5), "mm"))
  for (i in seq_len(nlevels(col_split)))
    decorate_annotation("A_score", slice = i,
      grid.lines(c(0, 1), unit(rep(thr["KMT2D"], 2), "native"), gp = gpar(lty = 2, col = "grey30", lwd = 0.7)))
  if (panelA$footnote)
    grid.text(sprintf("Controls: %d of %d held-out controls, randomly selected (seed %d).\nOther groups: all samples not used for CpG selection. Dashed line: threshold.",
                      length(ctrl_pick), length(ctrl_pool), panelA$seed),
              x = unit(2, "mm"), y = unit(1.5, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 5.5, col = "grey35"))
}
save_both_devices("heatmap_panelA_KMT2D", panelA$width, panelA$height, renderA)
