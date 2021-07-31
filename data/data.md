# The input data 

1. 在 `Linux Shell` 中，运行以下代码，下载 `input data`，必要时解压
2. 在浏览器中，输入以下 `url` ，下载 `input data`，并存储在 `data` 文件夹，必要时解压

```sh
mkdir data
```

```sh
# for: CellChat-vignette
# combined data from two biological conditions: normal and diseases
file="data_humanSkin_CellChat.rda"
url="https://ndownloader.figshare.com/files/25950872"
cd data
wget $url
cd ../
```

```sh
# for: Comparison_analysis_of_multiple_datasets
file1="cellchat_humanSkin_NL.rds"
url1="https://ndownloader.figshare.com/files/25954199"
file2="cellchat_humanSkin_LS.rds"
url2="https://ndownloader.figshare.com/files/25956518"
cd data
wget $url1
wget $url2
cd ../
```

```sh
# for: Comparison_analysis_of_multiple_datasets_with_different_cellular_compositions
file1="cellchat_embryonic_E13.rds"
url1="https://ndownloader.figshare.com/files/25957094"
file2="cellchat_embryonic_E14.rds"
url2="https://ndownloader.figshare.com/files/25957634"
cd data
wget $url1
wget $url2
cd ../
```
