# Analysis of mitochondrial genome
## Phylogenetic tree

We generated a (pseudo)-alignment of mitochondrial SNPs for all samples in the Europe+_recent dataset, including 5 rye owdery mildew genomes to be used as outgroup.

```
while read p; do
    bcftools consensus -a '-' --missing '-' --haplotype I -s ${p} -f ~/projects/project_data_prep/data/GCA_900519115.1_2022_bgt_ref_mating_type_MT880591.1.fa ~/projects/project_data_prep/data/2022+before2022+202
3+ncsu_covg15_recoded_snps_all_filtered_no_asterisk_MT880591.1.vcf.gz| sed "s/MT880591.1/${p}/g" > ${p}_mito_consensus.fa
done < list_Europe_recent_Bgs

cat *_mito_consensus.fa > Europe_recent_Bgs_mito.fa
rm *_mito_consensus.fa
```
We excluded all sites with more than 10% missing data with the script `variable_fasta.py" available from the repository [bion_tools](https://github.com/fmenardo/bion_tools).
```
python variable_fasta.py -m 0.9 Europe_recent_Bgs_mito.fa > Europe_recent_Bgs_mito_0.1_miss.fa
```
We inferred the phylogenetic tree with raxml-ng with 100 bootstrap replicate to calculate support.
```
raxml-ng --all --msa Europe_recent_Bgs_mito_0.1_miss.fa --model GTR+G+ASC_LEWIS --prefix Europe_recent_Bgs --threads 4 --workers 4 --bs-trees 100
```
The resulting tree was plotted with `plot_tree.R`.
