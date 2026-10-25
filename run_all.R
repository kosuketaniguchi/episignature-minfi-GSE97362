## =====================================================================
## run_all.R ― 01〜05 を順に実行する
## 初回は先に R/00_setup_packages.R と R/00_download_data.R を実行しておく。
## 各スクリプトは単独でも実行できる（すべて config.R を読む）。
## =====================================================================
steps <- c("R/01_gse97362_samplesheet.R",   # サンプルシートの作成
           "R/02_minfi_preprocess_qc.R",    # 前処理と QC（IDAT → β・M 行列）
           "R/03_build_episignatures.R",    # KMT2D・CHD7 シグネチャーの作成と検証
           "R/04_score_new_sample.R",       # 新しい検体の照合（デモ）
           "R/05_heatmap_episignatures.R")  # ヒートマップ
for (f in steps) {
  message("\n==== ", f, " (", format(Sys.time(), "%H:%M:%S"), ")")
  source(f, local = new.env())
}
message("\nすべて完了")
