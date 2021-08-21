# fastANCOM
A fast method for analysis of composition of microbiomes

## Install and load fastANCOM

```R
devtools::install_github('ZRChao/fastANCOM')
library(fastANCOM)
help(fastANCOM)
```

## News

- Version 0.0.1, 2021.7.20 Package first build 
- Version 0.0.2, 2021.8.01 We implement a fast version for the original ANCOM using the Mann-Whiteny U statistics, ```ANCOM(...,type='t2')```
- Version 0.0.3, 2021.8.10 Data pre-process step was implement, including detect structure zero, remove outlier, delete rare OTU or sample with low sequencing depth
- Version 0.0.4, 2021.8.21 Global and pairwise comparison among groups was designed
