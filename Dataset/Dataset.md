# Dataset
We uses the dataset <ins>Europe+</ins> from [Jigisha et al. 2024](https://doi.org/10.1101/2024.10.24.619980 ) [link](https://github.com/fmenardo/Bgt_popgen_Europe_2024/tree/Bgt_ms/Datasets). This datatset comprises 415 Bgt isolates that were sampled between 25 and 60 degree latitude North, and between -9 and 60 degree longitude.
For some analysis we used a subset of this dataset including only samples collected in 2015 or more recently (<ins>Europe+_recent</ins>). Finally we also defined a <ins>temporal</ins> dataset to investigate changes over multiple decades (see article for details). 

We use the vcf file produced by Jigisha et al. 2024 as a starting point for most of the analyses. The vcf can be downloaded [here](https://doi.org/10.5281/zenodo.13903934). For details about the WGS pipeline see the original repository [link](https://github.com/fmenardo/Bgt_popgen_Europe_2024/blob/Bgt_ms/WGS_pipeline/WGS_pipeline.md).

The list of all samples with metadata (including information about which sample was included in which dataset) can be found here [link](Minadakis_et_al_2024_Supplementary_Data_S1_rev.csv). Be aware there is a legend is at the bottom of the file, in case you want to import it e.g. in R you should remove it.
