# Manual of fastANCOM

A quick start for fastANCOM. fastANCOM is a simple and easy to operate R package, it aims to perform analysis of compositional microbiome data (ANCOM) at a fast speed.

# A toy example

Before general simulation comparison, we first present some results of a toy data generated from Poission to show how fastANCOM works. 
We called the results calculated from the log-ratio-transformed data as the **Joint model** and the results constructed based on the results from the log-transformed data as the **Marginal model**.
We shown that the consistence between the results calculated based on the marginal model and the results calculated based on the joint model directly.
We also shown that fastANCOM could provide a better effect size estimation than others.

# The simulation supplementary for fastANCOM
Simulated data for differential abundance (DA) analysis were generated from a Poisson- gamma mixture model. We used the code provided at https://github.com/FrederickHuangLin/ANCOM-BC-Code-Archive.
We evaluated the performance, in terms of power and false discovery rate (FDR), of fastANCOM with two zero-replacement strategies (Figures S1-S4). 
We then compared this version of fastANCOM with other DA testing methods, including edgeR, metagenomeSeq, ALDEx2, DESeq2, ANCOM, and ANCOM-BC. 

# Real data application

We applied fastANCOM to a human gut microbiome dataset from a study that investigated how fecal microbiomes differ across age and geography (Yatsunenko and others 2012. Human gut microbiome viewed across age and geography. Nature 486, 222–227). This dataset, which was also analyzed using ANCOM-BC, is available at https://github.com/FrederickHuangLin/ANCOM-BC-Code-Archive. It consists of bacterial V4 16S rRNA data from healthy children and adults from Amazonas of Venezuela, rural Malawi, and US metropolitan areas. For illustration, we carried out DA testing at the phylum level, separately for infants aged 0-2 years and adults aged 18-40 years, between each pair of populations (Malawi versus Venezuela, Malawi versus USA, and Venezuela versus USA).

