# fastANCOM
A fast method for analysis of composition of microbiomes

## Install and load fastANCOM

```R
devtools::install_github('ZRChao/fastANCOM')
library(fastANCOM)
help(fastANCOM)
```

## News

- 2021.7.20 Package first build
- 2021.8.01 We implement a fast version for the original ANCOM using the Mann-Whiteny U statistics, ```ANCOM(...,type='t2')```
