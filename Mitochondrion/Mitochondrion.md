# Analysis of mitochondrial genome
## Phylogenetic tree

We generated a (pseudo)-alignment of mitochondrial SNPs for all samples in the _Europe+_recent_ dataset, including 5 rye powdery mildew genomes to be used as outgroup.

```
while read p; do
    bcftools consensus -a '-' --missing '-' --haplotype I -s ${p} -f ~/projects/project_data_prep/data/GCA_900519115.1_2022_bgt_ref_mating_type_MT880591.1.fa ~/projects/project_data_prep/data/2022+before2022+202
3+ncsu_covg15_recoded_snps_all_filtered_no_asterisk_MT880591.1.vcf.gz| sed "s/MT880591.1/${p}/g" > ${p}_mito_consensus.fa
done < list_Europe_recent_Bgs

cat *_mito_consensus.fa > Europe_recent_Bgs_mito.fa
rm *_mito_consensus.fa
```
We excluded all sites with more than 10% missing data with the script `variable_fasta.py` available from the repository [bion_tools](https://github.com/fmenardo/bion_tools).
```
python variable_fasta.py -m 0.9 Europe_recent_Bgs_mito.fa > Europe_recent_Bgs_mito_0.1_miss.fa
```
We inferred the phylogenetic tree with raxml-ng with 100 bootstrap replicate to calculate support.
```
raxml-ng --all --msa Europe_recent_Bgs_mito_0.1_miss.fa --model GTR+G+ASC_LEWIS --prefix Europe_recent_Bgs --threads 4 --workers 4 --bs-trees 100
```
The resulting tree was plotted with `plot_tree.R`.

## Linkage disequilibrium
We calculated linkage disequilibrium as r2 between all pairs of variants within windows of 10kb over the mitochondrial and nuclear (11 chromosomes) genomes for the _Europe+_recent_ dataset using the script `launch_LD_decay_ext_eur_recent.sh`. To visualise LD decay with distance, we plotted the average r2 within distance classes of 100 bp against distance using the script `LD_decay_genome_mito.R`.

## Frequency of heterozygous calls
We obtained the number of high confidence calls (>90% read support) and the number of heterozygous calls (<90% read support) across the 11 chromosomes and the mitochondrial genome for all samples in the _Europe+_recent_ dataset. We used the output of the script [recode_multivcf_after_gatk_2024_new_clean.py](https://github.com/fmenardo/Bgt_popgen_Europe_2024/blob/Bgt_ms/WGS_pipeline/recode_multivcf_after_gatk_2024_new_clean.py) and plotted the distributions using the script `plot_het_comparsion.R`. 
