# episignature-minfi-GSE97362

DNA メチル化アレイ（Illumina 450K）の前処理から、エピシグネチャーによる照合までを R/Bioconductor の [minfi](https://bioconductor.org/packages/minfi/) で一通り行うコードです。
公開データ [GSE97362](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97362)（Kabuki 症候群・CHARGE 症候群、Butcher et al. 2017）を例に、KMT2D と CHD7 の 2 本のシグネチャーを作り、手元の 1 検体がどちらに当てはまるかを照合するところまでを再現します。

> 教育・研究用のコードです。臨床判断には、検証された検査を用いてください。

## 解析の流れ

| スクリプト | 内容 | 主な出力（`results/` の下） |
|---|---|---|
| `R/00_setup_packages.R` | 必要なパッケージのインストール（初回のみ） | ― |
| `R/00_download_data.R` | GSE97362 の IDAT・series matrix と、Zhou らのプローブ注釈を取得（初回のみ） | `data/GSE97362/` |
| `R/01_gse97362_samplesheet.R` | サンプルシートの作成（群・性別・年齢、発見／検証の区分） | `SampleSheet_minfi.csv` |
| `R/02_minfi_preprocess_qc.R` | 前処理と QC：IDAT から β 値・M 値の行列を作る | `minfi/` |
| `R/03_build_episignatures.R` | KMT2D・CHD7 シグネチャーの作成と検証 | `minfi_episignature/` |
| `R/04_score_new_sample.R` | 新しい検体の照合（GSE97362 の 2 検体を新患に見立てたデモ） | `minfi_new_sample/` |
| `R/05_heatmap_episignatures.R` | ヒートマップ（全体図、シグネチャーごとの教師なしクラスタリング、図用の小さな版） | `minfi_episignature/heatmap_*` |

## 使い方

```r
# 1. リポジトリを取得し、episignature-minfi-GSE97362.Rproj を RStudio で開く
#    （作業ディレクトリがリポジトリの最上位になる。コマンドラインなら最上位で Rscript を実行）

# 2. 初回だけ
source("R/00_setup_packages.R")
source("R/00_download_data.R")   # 数 GB。データを別の場所に置く場合は先に config.R を編集

# 3. 解析（01〜05 を順に実行）
source("run_all.R")
```

- パスの設定は `config.R` にまとめてあります。データや結果を別の場所に置くときは、`config.local.R` を作って `DATA_DIR`・`RESULTS_DIR` を書いてください（GitHub には上がりません）。環境変数 `GSE97362_DIR`・`RESULTS_DIR` でも指定できます。
- メモリは 8 GB 以上を推奨します。前処理と QC（02）は、235 検体でノート PC 10〜20 分程度です。
- 開発環境：macOS、R 4.6。

## 主な設計

**前処理と QC（02）**
- 検体の除外：検出 p < 0.01 のプローブが 95% 未満、全体強度（M・U の中央値の平均、log2）が 10.5 未満、バイサルファイト変換率が 80% 未満
- 身元の確認：推定性別と記録の不一致、SNP プローブによる混入の疑い（β 値の、遺伝型の山 0 / 0.5 / 1 からの平均距離が、コホートの中央値 + 5 MAD を超える）は解析から外し、QC 表に理由を残す
- 正規化：noob（背景と色素の補正）。検体ごとに独立に計算するので、新しい検体を参照コホートの再処理なしで照合できる
- プローブの除外：5% を超える検体で検出不良、Zhou らの推奨マスク（MASK_general）、CpG・単一塩基伸長位置の SNP（dbSNP）、交差反応、性染色体
- 細胞組成：IDOL（FlowSorted.Blood.EPIC）を採用。Houseman 法と EpiDISH は比較用

**シグネチャーの作成（03）**
- CpG の選抜には、GEO が定義する発見コホート（KMT2D LOF 例と対照、CHD7 LOF 例と対照）だけを使う
- M 値で limma（BH 調整 p < 0.05、|Δβ| ≥ 0.05）→ |t| の大きい順に、相関 0.85 を超える CpG を間引いて最大 200 個
- 共変量：KMT2D は年齢・性別のみ、CHD7 は細胞組成（5 分画）も追加。細胞組成を入れなかった場合に、対照の中でスコアとリンパ球割合が相関したかどうかで決めた（`check_cell_confounding_by_version.csv`）
- スコア：各 CpG の、発見コホートの対照からのずれ（z 値）を期待される方向に符号付けして平均。線形 SVM の確率も出す
- 陽性の閾値：選抜に使っていない対照のスコアの平均 + 3SD
- 選抜に使った検体で描いた図は、偶然の差でもきれいに分かれるので検証にならない（`heatmap_discovery_circular.pdf` で実演）

## 注意と限界

- 450K のデータで作ったシグネチャーです。EPIC v1・v2 の検体では一部の CpG が欠けます（04 は使えた CpG の数を表示します）。
- 性染色体は除外しています。X 連鎖の変化を見るときは、除外せずに男女別に解析してください。
- 参照は GSE97362 の対照（北米の小児コホート）です。別の集団（日本人など）の検体を照合するときは、集団差とバッチ差に注意してください。自施設の対照を同じチップで測り、閾値を超えないことを確かめるのが確実です。
- 図には日本語を含みます。macOS では quartz、それ以外では `pdf(family = "Japan1")` で出力します。

## 引用

このコードを使った場合は `CITATION.cff` の情報と、以下の原著を引用してください。

- データ：Butcher DT, et al. CHARGE and Kabuki syndromes: gene-specific DNA methylation signatures identify epigenetic mechanisms linking these clinically overlapping conditions. *Am J Hum Genet* 100:773–788, 2017
- minfi：Aryee MJ, et al. *Bioinformatics* 30:1363–1369, 2014
- noob：Triche TJ Jr, et al. *Nucleic Acids Res* 41:e90, 2013
- プローブ注釈：Zhou W, et al. *Nucleic Acids Res* 45:e22, 2017
- IDOL：Salas LA, et al. *Genome Biol* 19:64, 2018
- Houseman 法：Houseman EA, et al. *BMC Bioinformatics* 13:86, 2012
- EpiDISH：Teschendorff AE, et al. *BMC Bioinformatics* 18:105, 2017
- wateRmelon：Pidsley R, et al. *BMC Genomics* 14:293, 2013
- limma：Ritchie ME, et al. *Nucleic Acids Res* 43:e47, 2015
- ComplexHeatmap：Gu Z, et al. *Bioinformatics* 32:2847–2849, 2016

## ライセンス

コードは MIT ライセンスです（`LICENSE`）。GSE97362 のデータは GEO の利用条件に従います。

---

## English summary

R/Bioconductor (minfi) code for Illumina 450K DNA methylation arrays, from raw IDAT files to episignature scoring, using the public dataset GSE97362 (Kabuki and CHARGE syndromes; Butcher et al., *AJHG* 2017). The pipeline performs sample QC (detection, intensity, bisulfite conversion, sex check, SNP-based identity/contamination), noob normalization, probe masking (Zhou et al. 2017), and IDOL cell-composition estimation (`02`); builds KMT2D and CHD7 episignatures from the GEO-defined discovery cohorts with limma and validates them on held-out samples (`03`); scores new samples against the saved signatures without reprocessing the reference cohort (`04`); and draws heatmaps including unsupervised hierarchical clustering (`05`). Edit `config.R` to set data and output locations, then run `run_all.R`. Supplementary code for a Japanese review article in *Jikken Igaku* (Yodosha, 2027). For research and education only; not a clinical test.
