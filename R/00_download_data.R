## =====================================================================
## 00_download_data.R ― 解析に使う公開データを取得する（初回だけ。数 GB、時間がかかる）
##   1. GSE97362 の series matrix（群・性別・年齢などのメタデータ）
##   2. GSE97362 の RAW（IDAT）
##   3. Zhou らの 450K プローブ注釈（推奨マスク MASK_general と、集団別の SNP マスク）
## 置き場所は config.R の DATA_DIR。すでにあるファイルは取り直さない。
## =====================================================================
source("config.R")
suppressPackageStartupMessages(library(GEOquery))
options(timeout = max(3600, getOption("timeout")))
dir.create(IDAT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ANNOT_DIR, recursive = TRUE, showWarnings = FALSE)

## 1. series matrix
if (!file.exists(SERIES_FILE)) {
  getGEO("GSE97362", destdir = DATA_DIR, GSEMatrix = TRUE, getGPL = FALSE)
  if (!file.exists(SERIES_FILE)) stop("series matrix を保存できませんでした: ", SERIES_FILE)
}

## 2. RAW（IDAT の tar）。GEO の IDAT は .idat.gz のことがあるが、minfi も 01 もそのまま読める
if (!length(list.files(IDAT_DIR, pattern = "_Grn\\.idat(\\.gz)?$"))) {
  getGEOSuppFiles("GSE97362", makeDirectory = FALSE, baseDir = DATA_DIR, filter_regex = "RAW\\.tar$")
  tar_file <- list.files(DATA_DIR, pattern = "RAW\\.tar$", full.names = TRUE)[1]
  untar(tar_file, exdir = IDAT_DIR)
  message("展開した IDAT: ", length(list.files(IDAT_DIR, pattern = "\\.idat(\\.gz)?$")), " 本（235 検体 × 2 = 470 本のはず）")
}

## 3. Zhou らのプローブ注釈（https://zwdzwd.github.io/InfiniumAnnotation）
##    MASK_general 列があるのは 2022-09 のアーカイブ版
base <- "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/"
get_file <- function(url, dest) if (!file.exists(dest)) download.file(url, dest, mode = "wb")
get_file(paste0(base, "archive/202209/HM450.hg19.manifest.tsv.gz"), file.path(ANNOT_DIR, "HM450.hg19.manifest.tsv.gz"))
get_file(paste0(base, "HM450.hg19.manifest.pop.tsv.gz"),           file.path(ANNOT_DIR, "HM450.hg19.manifest.pop.tsv.gz"))
message("完了: ", normalizePath(DATA_DIR))
