# fastANCOM
A fast method for analysis of compositions of microbiomes. 

## Install and load fastANCOM

```R
devtools::install_github('ZRChao/fastANCOM')
library(fastANCOM)
help(fastANCOM)
```
The details of the manual could be found at https://rpubs.com/RChao/820127. The codes for simulation and real data application could be found at the fold ```../Supp```


## News

- Version 0.0.1, 2021.7.20 Package first build 
- Version 0.0.2, 2021.8.01 We implement a fast version for the original ANCOM using the Mann-Whiteny U statistics, ```ANCOM(...,type='t2')```
- Version 0.0.3, 2021.8.10 Data pre-process step was implement, including detect structure zero, remove outlier, delete rare OTU or sample with low sequencing depth
- Version 0.0.4, 2021.8.21 Global and pairwise comparison among groups was designed

Any suggestions or problem, please contact Chao Zhou（Supdream8@sjtu.edu.cn) .

*Citation*：
- Mandal, Siddhartha, et al. "Analysis of composition of microbiomes: a novel method for studying microbial composition." Microbial ecology in health and disease 26.1 (2015): 27663.
- Zhou, C., Wang, H., Zhao, H., & Wang, T. (2022). fastANCOM: a fast method for analysis of compositions of microbiomes. Bioinformatics, 38(7), 2039-2041.
