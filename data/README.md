# data/

解析に使う公開データの置き場所です。データそのものはリポジトリに含めません（`.gitignore`）。
`R/00_download_data.R` を実行すると、次の構成で自動的に取得されます。別の場所に置く場合は `config.R` の `DATA_DIR` を変えてください。

```
data/GSE97362/
├── GSE97362_series_matrix.txt.gz        # メタデータ（群・性別・年齢・発見／検証の区分）
├── idat/
│   ├── GSM2562699_9376525031_R05C01_Grn.idat
│   ├── GSM2562699_9376525031_R05C01_Red.idat
│   └── ...                              # 235 検体 × 2 = 470 本（.idat.gz のままでもよい）
└── annotation/
    ├── HM450.hg19.manifest.tsv.gz       # Zhou らの推奨マスク（MASK_general）
    └── HM450.hg19.manifest.pop.tsv.gz   # 集団別の SNP マスク（日本人検体では EAS 列）
```

- GSE97362：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97362
- プローブ注釈：https://zwdzwd.github.io/InfiniumAnnotation
