library(ggplot2)
library(dplyr)
library(stringr)

setwd("~/projects/nikos/LD_decay/2_output/10kb/")

# Europe+_recent MITOCHONDRIA

# import the data
mtld <- read.table("mt_LD_plink2_ext_eur_recent_10kb.vcor", sep="", header=F)

# check column names
colnames(mtld) <- c("CHROM_A", "BP_A", "ID_A", "CHROM_B", "BP_B", "ID_B", "R2")

mtld_summary_sort <- mtld %>%
  mutate(dist = BP_B - BP_A) %>%
  select(dist, R2) %>%
  arrange(dist)

mtld_summary_sort$dists <- cut(mtld_summary_sort$dist,
                             breaks=seq(from=min(mtld_summary_sort$dist)-1,
                                        to=max(mtld_summary_sort$dist)+1,
                                        by=100)) # by=10 makes it more detailed

mtdfr1 <- mtld_summary_sort %>% group_by(dists) %>% summarise(mean=mean(R2), median=median(R2))

mtdfr1 <- mtdfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))


# Europe+_recent 11 chromosomes

# import the data
ld <- read.table("LD_plink2_ext_eur_recent_10kb.vcor", sep="", header=F)

# check column names
colnames(ld) <- c("CHROM_A", "BP_A", "ID_A", "CHROM_B", "BP_B", "ID_B", "R2")

ld_summary_sort <- ld %>%
  mutate(dist = BP_B - BP_A) %>%
  select(dist, R2) %>%
  arrange(dist)

ld_summary_sort$dists <- cut(ld_summary_sort$dist,
                             breaks=seq(from=min(ld_summary_sort$dist)-1,
                                        to=max(ld_summary_sort$dist)+1,
                                        by=100)) # by=10 makes it more detailed

dfr1 <- ld_summary_sort %>% group_by(dists) %>% summarise(mean=mean(R2), median=median(R2))

dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(dists,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))


LD_decay <- ggplot() +
  geom_line(data=mtdfr1, aes(x=start, y=mean, colour="Mitochondrial"), linewidth=0.6, alpha=0.8)+
  geom_line(data=dfr1, aes(x=start, y=mean, colour="Nuclear"), linewidth=0.6, alpha=0.8)+
  labs(x="Distance (Kb)", y=expression(LD~(r^{2})))+
  scale_colour_manual(name = 'Genome', 
                      values = c('Mitochondrial'='red', 'Nuclear'='black')) +
  scale_x_continuous(breaks=c(0,2*10^3,4*10^3,6*10^3,8*10^3,10*10^3), labels=c("0","2","4","6","8","10"))+
  theme_bw()

pdf(paste("plots/mito_nuc_ext_eur_recent_10kb_by100.pdf"), width = 8, height = 6)
print(LD_decay)
dev.off()

png(paste("plots/mito_nuc_ext_eur_recent_10kb_by100.png"), width = 8, height = 6, units = "in", res = 150)
print(LD_decay)
dev.off()