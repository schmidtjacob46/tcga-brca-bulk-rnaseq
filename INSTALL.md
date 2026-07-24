# Software Requirements

## Python

Python 3.11 or newer

Install Python dependencies:

``` bash
pip install -r requirements.txt
```

## R

This project was developed using **R 4.2+** and requires the following
packages:

-   DESeq2
-   tidyverse
-   ggrepel

Install them with:

``` r
install.packages(c("tidyverse", "ggrepel"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```
