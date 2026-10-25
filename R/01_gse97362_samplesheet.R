## =====================================================================
## 01_gse97362_samplesheet.R  （オフライン版：ダウンロードはしない）
##
## 前提：config.R の DATA_DIR に（00_download_data.R で取得できる）
##   idat/GSM*_<SentrixID>_<R0xC0x>_{Grn,Red}.idat(.gz)   （235 検体 × 2 = 470 本。.gz のままでもよい）
##   GSE97362_series_matrix.txt.gz                  （メタデータはこのヘッダ部だけ使う）
## が既にある。GEOquery は使わず、series matrix の "!Sample_" 行だけを読む（数秒）。
##
## 出力：RESULTS_DIR/SampleSheet_minfi.csv  ―― minfi::read.metharray.exp() 用（Basename 列）
##   group : Control / KS1_KMT2D_LOF / KS1_KMT2D_P_LP / KS1_KMT2D_VUS / KS2_KDM6A /
##           CHARGE_CHD7_LOF / CHARGE_CHD7_P_LP / CHARGE_CHD7_VUS
##   set   : discovery_KMT2D（KMT2D LOF 11 例 + その対照 11 例）/ discovery_CHD7 / validation
##           ※ GEO の "sample type" に発見・検証の区分がそのまま入っている
## =====================================================================
suppressPackageStartupMessages(library(data.table))

source("config.R")
idat_dir <- IDAT_DIR
sm_file  <- SERIES_FILE
out_file <- SAMPLESHEET
if (!file.exists(sm_file)) stop("series matrix が見つかりません: ", sm_file, "（00_download_data.R を実行するか config.R の DATA_DIR を確認）")

## ---- 1. series matrix のヘッダ（!Sample_ 行）だけ読む ---------------------
con <- gzfile(sm_file, "rt")
hdr <- character(0)
while (length(l <- readLines(con, n = 1)) > 0) {
  if (startsWith(l, "!series_matrix_table_begin")) break
  if (startsWith(l, "!Sample_")) hdr <- c(hdr, l)
}
close(con)
split_row <- function(l) { x <- strsplit(l, "\t", fixed = TRUE)[[1]]; list(key = x[1], val = gsub('^"|"$', "", x[-1])) }
rows <- lapply(hdr, split_row)
get1 <- function(key) rows[[which(sapply(rows, `[[`, "key") == key)[1]]]$val
gsm   <- get1("!Sample_geo_accession")
title <- get1("!Sample_title")           # 例 "CHD7-1 whole blood CHD7 LOF genomic DNA"
desc  <- get1("!Sample_description")     # 例 "CHD7-1"（個体ID。Sample_Name に使う）
n <- length(gsm)

## characteristics（gender / age (years) / sample type / disease state / tissue）を列に展開
chars <- rows[sapply(rows, `[[`, "key") == "!Sample_characteristics_ch1"]
char_dt <- rbindlist(lapply(chars, function(r) {
  kv <- tstrsplit(r$val, ": ", fixed = TRUE, keep = 1:2)
  data.table(i = seq_len(n), field = kv[[1]], val = kv[[2]])   # 列名 "key" は data.table() の予約引数なので不可
}))
char_wide <- dcast(char_dt, i ~ field, value.var = "val")

## ---- 2. ラベルの整形 -------------------------------------------------
variant_class <- sub("^\\S+ whole blood (.+) genomic DNA$", "\\1", title)   # "Control", "KMT2D LOF", "CHD7 VUS", ...
to_group <- function(vc) {
  if (vc == "Control") return("Control")
  g <- sub(" .*$", "", vc); cls <- sub("^\\S+ ", "", vc)
  cls <- switch(cls, "LOF" = "LOF", "VUS" = "VUS", "pathogenic" = "P_LP", "likely pathogenic" = "P_LP", cls)
  if (g == "KDM6A") return("KS2_KDM6A")
  paste0(c(KMT2D = "KS1_KMT2D", CHD7 = "CHARGE_CHD7")[[g]], "_", cls)
}
to_set <- function(st) fcase(grepl("KMT2D LOF discovery", st), "discovery_KMT2D",
                             grepl("CHD7 LOF discovery",  st), "discovery_CHD7",
                             default = "validation")
parse_age <- function(a) {              # "NA" → NA、"<0.1" → 0.05（生後1か月未満）
  a <- trimws(a); out <- suppressWarnings(as.numeric(a))
  lt <- grepl("^<", a); out[lt] <- as.numeric(sub("^<", "", a[lt])) / 2
  out
}

ss <- data.table(
  Sample_Name   = desc,
  gsm           = gsm,
  group         = unname(vapply(variant_class, to_group, character(1))),
  set           = to_set(char_wide[["sample type"]]),
  disease_state = char_wide[["disease state"]],
  sample_type   = char_wide[["sample type"]],
  variant_class = variant_class,
  gene          = sub("-.*$", "", desc),
  Sex           = c(male = "M", female = "F")[char_wide[["gender"]]],
  Age           = parse_age(char_wide[["age (years)"]]),
  Age_raw       = char_wide[["age (years)"]]
)
stopifnot(!anyDuplicated(ss$Sample_Name))

## ---- 3. IDAT の Basename（_Grn/_Red.idat を除いたファイル名の幹）を GSM に紐づける ----
grn <- list.files(idat_dir, pattern = "^GSM\\d+_\\d+_R\\d+C\\d+_Grn\\.idat(\\.gz)?$")
if (!length(grn)) stop("IDAT が見つかりません: ", idat_dir)
bn <- data.table(gsm = sub("_.*$", "", grn),
                 Sentrix_ID = sub("^GSM\\d+_(\\d+)_R\\d+C\\d+_Grn\\.idat(\\.gz)?$", "\\1", grn),
                 Sentrix_Position = sub("^GSM\\d+_\\d+_(R\\d+C\\d+)_Grn\\.idat(\\.gz)?$", "\\1", grn),
                 Basename = sub("_Grn\\.idat(\\.gz)?$", "", grn))   # ファイル名の幹だけを保存（場所は 02・04 で IDAT_DIR から組み立てる）
red_ok <- file.exists(file.path(idat_dir, paste0(bn$Basename, "_Red.idat"))) | file.exists(file.path(idat_dir, paste0(bn$Basename, "_Red.idat.gz")))
stopifnot(all(red_ok))   # Red が揃っているか
ss <- merge(ss, bn, by = "gsm", all.x = TRUE, sort = FALSE)
if (anyNA(ss$Basename)) stop("IDAT が見つからない GSM: ", paste(ss$gsm[is.na(ss$Basename)], collapse = ", "))

setcolorder(ss, c("Sample_Name", "gsm", "group", "set", "disease_state", "sample_type", "variant_class", "gene",
                  "Sex", "Age", "Age_raw", "Sentrix_ID", "Sentrix_Position", "Basename"))
fwrite(ss, out_file)

## ---- 4. 確認表示（期待値：235 検体、discovery_KMT2D = LOF 11 + Control 11、chip 38 枚）----
cat("書き出し:", out_file, "\n検体数:", nrow(ss), " chip 数:", uniqueN(ss$Sentrix_ID), "\n")
print(dcast(ss[, .N, by = .(group, set)], group ~ set, value.var = "N", fill = 0))
cat("Age 欠損/不定:", paste(ss[is.na(Age) | grepl("^<", Age_raw), paste0(Sample_Name, "(", Age_raw, ")")], collapse = ", "), "\n")
print(ss[set == "discovery_KMT2D", .(n = .N, age_median = median(Age), age_min = min(Age), age_max = max(Age)), by = group])
