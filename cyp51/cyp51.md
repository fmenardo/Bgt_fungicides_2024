# Cyp51 coverage and alignmnent

We estimated the number of copies of cyp51 by calculating the average coverage from start to stop codon and dividing it by the genome-wide coverage.
We used samtools v1.17, for example, for one isolate:
```
samtools coverage -r LR026991.1_chr8:5728014-5729706 ${filename}
```
The genome-wide covergae was obtained from (LINK). The coverages and the ratio are availbles in ` cyp51_cnv.csv`

Moreover we obtained a consensus sequence in which we reported all mutations, also heterozygous ones. Therefore in this alignment a mutation is present even if it is present in only one copy of cyp51.

```

bcftools view  -o cyp51_het.vcf.gz -O z -a -R ../Fungicides_targets/coord_cyp51 -S tritici_extended_europe_2022+before2022+2023+ncsu -f GCA_900519115.1_2022_bgt_ref_mating_type_LR026991.1_chr8.fa 2022+before2022+2023+ncsu_covg15_LR026991.1_chr8.vcf.gz
bcftools index cyp51_het.vcf.gz

python make_fasta_gene_from_vcf.py -vcf cyp51_het.vcf.gz -o Bgt_Eur+_cyp51_het -rc

```


