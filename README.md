# Statistical-analyses-optimal-sampling-size-in-Mediterranean-forest-ecosystems
Simulated data and R code to perform the statistical analyses to identify the minimum number of pooled samples needed to obtain a reliable description of fungal communities in terms of diversity and composition in three different Mediterranean forests' stands. This is the code used in Adamo et al. to perform the statistical analyses 
that can be used in similar works. The repository was published in Zenodo (doi:10.5281/zenodo.4434407).

By using high throughput-DNA sequencing data (Illumina MiSeq), we identified the minimum number of pooled samples needed to obtain a reliable description of fungal communities
in terms of diversity and composition in three different Mediterranean forests' stands. Twenty soil samples were randomly taken and kept separated in each of the three plots per forest type. After keeping one random sample which was not pooled, we obtained 5 composite samples consisting of pools of 3, 6, 10, 15 and 20 samples. We then sequenced the ITS2 using Illumina MiSeq.
We further tested:
1. The effect of sample pooling on fungal diversity, the iNEXT function was used to build rarefactions curves pooling together the individual samples. 

2. The variance of Bray-Curtis matrix between the number of sample pools for each forest type was compared using the betadisper function which is analogue to a Levene's test.

3. Beta-diversity patterns and whether the core of most abundant fungal species is maintained between sites, we evaluated for each pool the species (or abundances-per-species) 
losses (B) and species gains (C) using the beta-indices (tbi function, Legendre, 2019, Ecology and Evolution). We used one-sample pool per each forest (sample 1) as a reference, and we compared pools with increasing number of samples (sample 3, 6, 10, 15 and 20) to identify species losses and gains.

The results of these analyses are included in : I. Adamo, Y. Pinuela, J.A. Bonet, C. Castaño, J. Martínez de Aragón, J. Parladé, J. Pera, J.G. Alday, 2021. Sampling forest soils to describe fungal diversity and composition. Which is the optimal sampling size in mediterranean pure and mixed pine oak forests?, Fungal Biology, https://doi.org/10.1016/j.funbio.2021.01.005

This project has received funding from European's Union H2020 research and innovation programme under Marie Slodowska-Curie grant agreement No 801586.
