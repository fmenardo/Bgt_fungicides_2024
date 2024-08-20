# cytb

We generated a SNP alignment with all mitochondrial SNPs with less than 10% of missing data. The script variable_fasta.py is available [here](https://github.com/fmenardo/bion_tools).

```
while read p; do
    bcftools consensus -a '-' --missing '-' --haplotype I -s ${p} -f ../../project_data_prep/data/GCA_900519115.1_2022_bgt_ref_mating_type_MT880591.1.fa ../../project_data_prep/data/2022+before2022+2023+ncsu_covg15_recoded_snps_all_filtered_no_asterisk_MT880591.1.vcf.gz| sed "s/MT880591.1/${p}/g" > ${p}_mito_consensus.fa
done < list_Bgt_Europe+_Bgs

cat *_mito_consensus.fa > mito_msa.fa
rm *_mito_consensus.fa

python variable_fasta.py -m 0.9 mito_msa.fa > mito_msa_10%miss_max.fa
```

