setwd("~/projects/project_data_prep/data/")

library(ggplot2)
library(patchwork, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.2/")

gen_len = 137402123
n_pos = list(5800,9600,2600,5800,8400,8035,7800,5200,10000,6000,2000)
wg_N = Reduce('+',n_pos)   ## 71235

wg_stats <- read.csv("../analysis/pl_stats/2022+before2022+2023+ncsu_covg15_chr1-11_recoding_stats_17042024.csv")
wg_stats$WG_pos_without_N <- wg_stats$WG_TotalNumPos - wg_N
wg_stats$WG_high_conf <- wg_stats$WG_pos_without_N - wg_stats$WG_NumMissingPos 
wg_stats$conf_het_ratio <-  wg_stats$WG_NumHetPos / wg_stats$WG_high_conf

ggplot(data=wg_stats)+geom_histogram(aes(conf_het_ratio),fill='plum3',alpha=0.4,color='plum3')+
  theme_minimal()+labs(x="Number of high confidence calls / Number of heterozygous positions")

ext_eur_rec <- read.csv("~/projects/vcf_project_tritici/tritici_recent_extended_europe_2022+2023+ncsu.args", header = FALSE)
wg_eur <- merge(wg_stats, ext_eur_rec, by.x = "Isolate", by.y = "V1")

ggplot(data=wg_eur)+geom_histogram(aes(conf_het_ratio),fill='plum3',alpha=0.4,color='plum3')+
  theme_minimal()+labs(x="Number of high confidence calls / Number of heterozygous positions")

mt_stats <- read.csv("2022+before2022+2023+ncsu_covg15_recoded_MT880591.1_stats.csv")
mt_stats$high_conf <- 109707 - mt_stats$NumMissingPos
mt_stats$conf_het_ratio <- mt_stats$NumHetPos / mt_stats$high_conf
mt_eur <- merge(mt_stats, ext_eur_rec, by.x = "Isolate", by.y = "V1")

confhetpos <- ggplot()+geom_histogram(data=wg_eur, aes(conf_het_ratio, fill='Nuclear',colour='Nuclear'),alpha=0.5,)+
  geom_histogram(data=mt_eur, aes(conf_het_ratio, fill='Mitochondrial', colour = 'Mitochondrial'), alpha = 0.5) +
  scale_fill_manual(name = 'Genome',
                    values = c('Nuclear'='black','Mitochondrial'='red'))+
  scale_colour_manual(name = 'Genome', 
                      values = c('Nuclear'='black','Mitochondrial'='red'))+theme_bw()+
  labs(x="Number of heterozygous sites / Number of high confidence sites",y="Count")+
  coord_cartesian(xlim=c(0,0.007))

pdf(paste("~/projects/project_bgt_mt_tree/analysis/NEW_hetpos_comp_blackred.pdf"), width = 8, height = 6)
print(confhetpos)
dev.off()
