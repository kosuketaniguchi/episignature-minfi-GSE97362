## =====================================================================
## 00_setup_packages.R ― 必要なパッケージを入れる（初回だけ）
## R >= 4.3 / Bioconductor >= 3.18 を想定（EPIC v2 の注釈パッケージは 3.18 以降）
## =====================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c(
  "minfi",                                         # 本体：IDAT の読込・noob・QC・性別推定・SNP プローブ
  "IlluminaHumanMethylation450kmanifest",          # GSE97362（450K）用のマニフェスト
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",  # 450K の注釈（hg19）
  "IlluminaHumanMethylationEPICmanifest",          # EPIC v1（自分の検体が EPIC のとき）
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICv2manifest",        # EPIC v2
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
  "wateRmelon",             # bscon()（バイサルファイト変換率）、BMIQ()（プローブの型の補正、任意）
  "FlowSorted.Blood.EPIC",  # estimateCellCounts2()（IDOL 参照による血球組成の推定）
  "FlowSorted.Blood.450k",  # Houseman 法による推定（IDOL との比較用）
  "EpiDISH",                # 細胞組成推定の代替（β 行列だけで動く）
  "limma",                  # M 値の群間比較
  "ComplexHeatmap",         # ヒートマップ（05）
  "GEOquery"                # データの取得（00_download_data.R）
)
BiocManager::install(setdiff(bioc_pkgs, rownames(installed.packages())), update = FALSE, ask = FALSE)

cran_pkgs <- c("data.table", "jsonlite", "e1071", "pheatmap", "ggplot2", "RColorBrewer", "circlize")
install.packages(setdiff(cran_pkgs, rownames(installed.packages())))

## 任意：交差反応プローブのリスト（Chen 2013 / Pidsley 2016 / McCartney 2016）
## remotes::install_github("markgene/maxprobes")
## 任意：SNP プローブによる同一性・混入の詳細な確認
## remotes::install_github("hhhh5/ewastools")
