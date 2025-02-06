#!/bin/bash

#SBATCH --time=0-24:00:00
#SBATCH --mem 10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --job-name=LD_decay_plink2_10kb



#####change the args and names


##### europe+_recent dataset - 11 chromosome

# gatk filter and subset (FILTERING FOR AC is not really needed because plink only measures correlations between actual variants, but we do it to reduce VCF size)
~/data/gatk-4.4.0.0/gatk SelectVariants -V ~/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
-O ~/projects/nikos/LD_decay/0_data/tritici_ext_eur_recent_biallelic_snps_filtered_miss0.1.vcf.gz \
-sn ~/projects/vcf_project_tritici/tritici_recent_extended_europe_2022+2023+ncsu.args \
--exclude-intervals MT880591.1 \
--exclude-intervals LR026995.1_Un \
--exclude-intervals Bgt_MAT_1_1_3 \
--max-nocall-fraction 0.1 \
--select "AC > 0 && AC < 368" 

# vcf to plink format
# whole dataset
~/data/plink2 --vcf ~/projects/nikos/LD_decay/0_data/tritici_ext_eur_recent_biallelic_snps_filtered_miss0.1.vcf.gz \
--make-bed --out ~/projects/nikos/LD_decay/0_data/tritici_ext_eur_recent_biallelic_snps_filtered_miss0.1

# run LD decay ld-window 10 kb
# whole dataset
~/data/plink2 --bfile ~/projects/nikos/LD_decay/0_data/tritici_ext_eur_recent_biallelic_snps_filtered_miss0.1 \
--r2-unphased --ld-window-r2 0 --ld-window-kb 10 --out ~/projects/nikos/LD_decay/2_output/10kb/LD_plink2_ext_eur_recent_10kb

###### mitochondria ######
#### europe+_recent dataset
~/data/gatk-4.4.0.0/gatk SelectVariants -V ~/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
-O ~/projects/nikos/LD_decay/0_data/mt_tritici_ext_eur_recent_biallelic_snps_filtered_miss0.1.vcf.gz \
-sn ~/projects/vcf_project_tritici/tritici_recent_extended_europe_2022+2023+ncsu.args \
--intervals MT880591.1 \
--max-nocall-fraction 0.1 \
--select "AC > 0 && AC < 368"

# vcf to plink format
# whole dataset
~/data/plink2 --vcf ~/projects/nikos/LD_decay/0_data/mt_tritici_ext_eur_recent_biallelic_snps_filtered_miss0.1.vcf.gz \
--make-bed --out ~/projects/nikos/LD_decay/0_data/mt_tritici_ext_eur_recent_biallelic_snps_filtered_miss0.1

# run LD decay ld-window 10 kb
# whole dataset
~/data/plink2 --bfile ~/projects/nikos/LD_decay/0_data/mt_tritici_ext_eur_recent_biallelic_snps_filtered_miss0.1 \
--r2-unphased --ld-window-r2 0 --ld-window-kb 10 --out ~/projects/nikos/LD_decay/2_output/10kb/mt_LD_plink2_ext_eur_recent_10kb
