# isoRelate analysis of cyp51 locus
We found previously that the cyp51 locus had an excess of Identity by descent (IBD) segments between pairs of isolates in many populations [link](https://github.com/fmenardo/Bgt_popgen_Europe_2024/blob/main/isoRelate/isoRelate.md).

To investigate this further we perfomrmed a new isoRelate analysis on this locus using all the isolates in the [Europe+_recent dataset](https://github.com/fmenardo/Bgt_popgen_Europe_2024/blob/main/Datasets/Datasets.md).
The scripts to run isoRelate and  the list of SNPs that can unambiguosly positions onto the genetic map `THUN12x96224_genetic_map_in_cM_+_phy_distance.ambiguous_intervals.list` are available [here](https://github.com/fmenardo/Bgt_popgen_Europe_2024/blob/main/isoRelate/isoRelate.md)

First we selected SNPs in a region of 3 mb around the gene cyp51

```
gatk SelectVariants \
     -R ~/projects/vcf_project_tritici/GCA_900519115.1_2022_bgt_ref_mating_type.fa \
     -V ~/projects/vcf_project_tritici/tritici_2022+before2022+2023+ncsu_ALL_biallelic_snps.vcf.gz \
     -O Europe+_recent_chr8_cyp51.vcf.gz \
     --max-nocall-fraction 0 \
     --exclude-intervals ~/projects/project_tritici_fabrizio/data/THUN12x96224_genetic_map_in_cM_+_phy_distance.ambiguous_intervals.list \
     --intervals LR026991.1_chr8:4000000-7000000 \
     --select "AF > 0.05 && AF < 0.95" \
     --sample-name ~/projects/vcf_project_tritici/tritici_recent_extended_europe_2022+2023+ncsu.args
```
We generated the input files for isoRelate

```
plink --allow-extra-chr --vcf Europe+_recent_chr8_cyp51.vcf.gz --recode 12 --double-id --out BgtE+r_chr8_cyp51 --threads 1

awk '$5="1" && $6="0"' BgtE+r_chr8_cyp51.ped >  BgtE+r_chr8_cyp51_mod.ped
sed -E 's/\S+_chr//g' BgtE+r_chr8_cyp51.map > BgtE+r_chr8_cyp51_mod.map


python add_cM_to_map.py -map BgtE+r_chr8_cyp51_mod.map -rec ~/projects/project_tritici_fabrizio/data/THUN12x96224_genetic_map_in_cM_+_phy_distance -o BgtE+r_chr8_cyp51_mod
```
We then inferred IBD segments with isoRelate

```
Rscript run_ibd_step1.R -o BgtE+r_cyp51 -p BgtE+r_chr8_cyp51_mod.ped -m BgtE+r_chr8_cyp51_mod_cM.map -c 5 -C 2
```
Finally we identified clusters, and plotted the graphs with `plot_clusters.R`, while the plot of the coverage of the tandem repeat was produced with `plot_coverage_dup.R`
