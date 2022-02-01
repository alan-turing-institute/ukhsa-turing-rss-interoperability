# Case study 2


## Getting started

To run the analysis requires installing the [R-INLA](https://www.r-inla.org) package:

```{R}
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
```

## Getting started

1. [prep_data.R](prep_data.R) script:
    - combines prevalence estimates with the census and IMD data to produce `prev_data.rds` for analysis
    - creates `data/W.adj` and `data/space_obj.Rdata`, which are used for handling spatial components in the analysis and plotting
2. [analysis.R](analysis.R) script:
    - fit both models and save the outputs to file
3. [plot.R](plot.R) script:
    - retrieve results 
    - create plots 


## Data overview

The data files include:
- [Prevalence estimates](data/logit_moments.csv) obtained using the [prevdebiasr](https://github.com/alan-turing-institute/prevdebiasr) package as described in case study 1
- Ethnicity data from the 2011 census
- Index of Multiple Deprivation by LTLA
- UK geography files (e.g., mapping LSOAs to LTLAs, LADs to Regions, etc.)