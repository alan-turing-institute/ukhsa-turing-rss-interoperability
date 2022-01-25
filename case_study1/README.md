# README for case study 1

Currently reproducibility pipeline for case_study1:

1. Clone repository from command line  
`git clone https://github.com/alan-turing-institute/ukhsa-turing-rss-interoperability.git`
3. Download [data_from_archive] (https://www.dropbox.com/s/xvfflubc2wqmpxt/data_from_archive.zip?dl=1) and unpack to 'case_study1/data_from_archive'
4. In an R session, install the renv R package if you don't have it already (e.g. via `install.packages("renv")`). Then, run the following
```
setwd('case_study1')
renv::activate()
renv::restore()
source('scripts/06_gather_IR_data_and_create_plots.R')
```


