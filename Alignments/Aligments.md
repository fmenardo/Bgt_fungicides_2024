# Alignments

First we create a list with all samples including 5 Bgs strains as outgroup
`
cat ../Dataset/tritici_extended_europe_2022+before2022+2023+ncsu.args ../Dataset/tritici_outgroups.args > tritici_extended_europe_2022+before2022+2023+ncsu+Bgs
`

For each target we extracted a vcf file with all positions covering the coding sequence with bcftools:

```
bcftools view  -o sdhB.vcf.gz -O z -a -R ../Fungicides_targets/coord_sdhB -S tritici_extended_europe_2022+before2022+2023+ncsu+Bgs -f GCA_900519115.1_2022_bgt_ref_mating_type_LR026992.1_chr9.fa 2022+before2022+2023+ncsu_covg15_recoded_LR026992.1_chr9.vcf.gz 
srun bcftools index sdhB.vcf.gz

bcftools view  -o sdhC.vcf.gz -O z -a -R ../Fungicides_targets/coord_sdhC -S tritici_extended_europe_2022+before2022+2023+ncsu+Bgs -f GCA_900519115.1_2022_bgt_ref_mating_type_LR026986.1_chr3.fa 2022+before2022+2023+ncsu_covg15_recoded_LR026986.1_chr3.vcf.gz
bcftools index sdhC.vcf.gz

bcftools view  -o sdhD.vcf.gz -O z -a -R ../Fungicides_targets/coord_sdhD -S tritici_extended_europe_2022+before2022+2023+ncsu+Bgs -f GCA_900519115.1_2022_bgt_ref_mating_type_LR026985.1_chr2.fa 2022+before2022+2023+ncsu_covg15_recoded_LR026985.1_chr2.vcf.gz
bcftools index sdhD.vcf.gz

bcftools view  -o cyp51.vcf.gz -O z -a -R ../Fungicides_targets/coord_cyp51 -S tritici_extended_europe_2022+before2022+2023+ncsu -f GCA_900519115.1_2022_bgt_ref_mating_type_LR026991.1_chr8.fa 2022+before2022+2023+ncsu_covg15_recoded_LR026991.1_chr8.vcf.gz
srun bcftools index cyp51.vcf.gz

bcftools view  -o Btub.vcf.gz -O z -a -R ../Fungicides_targets/coord_Btub -S tritici_extended_europe_2022+before2022+2023+ncsu+Bgs -f GCA_900519115.1_2022_bgt_ref_mating_type_LR026993.1_chr10.fa 2022+before2022+2023+ncsu_covg15_recoded_LR026993.1_chr10.vcf.gz
bcftools index Btub.vcf.gz

bcftools view  -o cytB.vcf.gz -O z -a -R ../Fungicides_targets/coord_cytB -S tritici_extended_europe_2022+before2022+2023+ncsu+Bgs -f GCA_900519115.1_2022_bgt_ref_mating_type_MT880591.1.fa 2022+before2022+2023+ncsu_covg15_recoded_MT880591.1.vcf.gz
bcftools index cytB.vcf.gz
```
We converts each vcf in a multiple sequence alignment fasta file, translate to proteins and convert using biopython:

```
python make_fasta_gene_from_vcf.py -vcf sdhB.vcf.gz -o Bgt_Eur+_sdhB -rc
python make_fasta_gene_from_vcf.py -vcf sdhC.vcf.gz -o Bgt_Eur+_sdhC -rc
python make_fasta_gene_from_vcf.py -vcf sdhD.vcf.gz -o Bgt_Eur+_sdhD
python make_fasta_gene_from_vcf.py -vcf cyp51.vcf.gz -o Bgt_Eur+_cyp51 -rc
python make_fasta_gene_from_vcf.py -vcf Btub.vcf.gz -o Bgt_Eur+_Btub
#python make_fasta_gene_from_vcf.py -vcf cytB.vcf.gz -o Bgt_Eur+_cytB -mt
```


