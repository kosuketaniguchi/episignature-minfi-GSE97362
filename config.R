## =====================================================================
## config.R ― パスの設定（自分の環境に合わせて変えるのはここだけ）
##
## すべてのスクリプトは、このファイルがある場所（リポジトリの最上位）を作業ディレクトリにして実行する。
## RStudio なら episignature-minfi-GSE97362.Rproj を開けばそうなる。
## 環境変数 GSE97362_DIR と RESULTS_DIR でも上書きできる。
##
## 自分の環境だけの設定は config.local.R に書く（GitHub には上がらない。.gitignore 済み）。例：
##   DATA_DIR    <- "~/data/GSE97362"   # データを別の場所に置いている
##   RESULTS_DIR <- DATA_DIR            # 結果もデータと同じ場所に書く
## =====================================================================
DATA_DIR    <- Sys.getenv("GSE97362_DIR", unset = file.path("data", "GSE97362"))  # IDAT と series matrix
RESULTS_DIR <- Sys.getenv("RESULTS_DIR",  unset = "results")                      # 解析結果
if (file.exists("config.local.R")) source("config.local.R")
DATA_DIR <- path.expand(DATA_DIR); RESULTS_DIR <- path.expand(RESULTS_DIR)

## ---- ここから下は通常は変えない ----
IDAT_DIR      <- file.path(DATA_DIR, "idat")                            # GSM*_{Grn,Red}.idat(.gz)
SERIES_FILE   <- file.path(DATA_DIR, "GSE97362_series_matrix.txt.gz")   # メタデータ（群・性別・年齢）
ANNOT_DIR     <- file.path(DATA_DIR, "annotation")                      # Zhou らのプローブ注釈（マスク）
SAMPLESHEET   <- file.path(RESULTS_DIR, "SampleSheet_minfi.csv")        # 01 の出力
PREPROC_DIR   <- file.path(RESULTS_DIR, "minfi")                        # 02 の出力（β・M 行列、QC）
SIG_DIR       <- file.path(RESULTS_DIR, "minfi_episignature")           # 03・05 の出力
NEWSAMPLE_DIR <- file.path(RESULTS_DIR, "minfi_new_sample")             # 04 の出力
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
