# Copy number variation

We estimated the number of copies of all target genes by calculating the average coverage from start to stop codon and dividing it by the genome-wide coverage.

We calculated coverage for the target genes and for the mitochondrion with samtools v1.17:
```
rm cov_cyp51_raw
for i in *.bam; do echo $i >> cov_cyp_51raw; samtools coverage -r LR026991.1_chr8:5728014-5729706 $i | cut -f 7 >> cov_cyp_51raw; done;

rm cov_Btub_raw
for i in *.bam; do echo $i >> cov_Btub_raw; samtools coverage -r LR026993.1_chr10:7766015-7767656 $i | cut -f 7 >> cov_Btub_raw; done;

rm cov_erg24_raw
for i in *.bam; do echo $i >> cov_erg24_raw; samtools coverage -r LR026990.1_chr7:6622512-6624033 $i | cut -f 7 >> cov_erg24_raw; done;

rm cov_erg2_raw
for i in .bam; do echo $i >> cov_erg2_raw; samtools coverage -r LR026988.1_chr5:18537698-18538493 $i | cut -f 7 >> cov_erg2_raw; done;

rm cov_cytB_raw
for i in .bam; do echo $i >> cov_cytB_raw; samtools coverage -r MT880591.1:97534-102064 $i | cut -f 7 >> cov_cytB_raw; done;

rm cov_mit_raw
for i in *.bam; do echo $i >> cov_mit_raw; samtools coverage -r MT880591.1 $i | cut -f 7 >> cov_mit_raw; done;

rm cov_sdhB_raw
for i in *.bam; do echo $i >> cov_sdhB_raw; samtools coverage -r LR026992.1_chr9:5629949-5630812 $i | cut -f 7 >> cov_sdhB_raw; done;

rm cov_sdhC_raw
for i in *.bam; do echo $i >> cov_sdhC_raw; samtools coverage -r LR026986.1_chr3:3466484-3467217 $i | cut -f 7 >> cov_sdhC_raw; done;

rm cov_sdhD_raw
for i in *.bam; do echo $i >> cov_sdhD_raw; samtools coverage -r LR026985.1_chr2:13068412-13069018 $i | cut -f 7 >> cov_sdhD_raw; done;
```
The genome-wide coverage was obtained from (LINK). 

We compiled all coverages and ratios with the R script `coverage_fungicides_targets.R`. The coverages and the ratio are availbles in `final_ds_table.csv`
