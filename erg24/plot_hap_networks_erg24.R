setwd("~/projects/project_fungicides/analysis/isoRelate_erg24/")
library(adegenet)
library(pegas)
library(seqinr)
library(dplyr)
library(tidyr)
library(ape)
library(seqinr)

#load("hap_net_xy.Rdata")  # pre cooked coordinates for hap net plotting, can be commented out, plot will be less pretty


seq<-read.fasta("~/projects/project_fungicides/analysis/Bgt_Eur+_erg24.fa",as.string=TRUE)
names<-c()
sequences <-list()
fasta_names <-list()
for (i in 1:length(seq)){
  if (grepl("n", seq[[i]][1])){   ## elimintae seq with N
    next}
  if (grepl("S-",attr(seq[i][[1]], "name"))){ ### eliminate Bgs
    next}
      sequences <- c(sequences, seq[[i]][1])
      fasta_names  <- c(fasta_names,attr(seq[i][[1]], "name"))
      names <- c(names,attr(seq[i][[1]], "name"))
}

write.fasta(sequences, fasta_names, "Bgt_Eur+_erg24_no_N.fa", as.string = TRUE)

ali<-fasta2DNAbin("Bgt_Eur+_erg24_no_N.fa")



names1<-read.table("~/projects/vcf_project_tritici/tritici_extended_europe_2022+before2022+2023+ncsu.args")
names1

names1 <- names1  %>% 
  filter(names1[,1] %in% names)

## read metainfo

meta<-read.csv("~/projects/vcf_project_tritici/2022+before2022+2023+ncsu_metadata+fs+admxK9_03052024.csv")
meta1<-data.frame(meta$Sample.Name,meta$Country,meta$fs_level_4,meta$Longitude,meta$Latitude,meta$Year.of.Collection)
colnames(meta1) <- c("Sample.Name","Country", "Population","Longitude","Latitude","Year")



df <- meta1  %>% 
  filter(Sample.Name %in% names1[,1])

hap_table<-read.csv("~/projects/nikos/fungicide/resistance/2_output/final_metadata_9_7_24.csv")
hap1 <- data.frame(hap_table$Sample.Name,hap_table$erg24_Y165F,hap_table$erg24_V295L,hap_table$erg24_F289H)

hap<-merge(df,hap1,by.x = "Sample.Name", by.y = "hap_table.Sample.Name")
colnames(hap) <- c("Sample.Name","Country", "Population","Longitude","Latitude","Year","Y165F","V295L","F289H")
hap <- hap %>%
  mutate(haplotype = paste(Y165F, V295L, F289H, sep = ""))






h <- haplotype(ali)
(P <- haploFreq(ali, fac = hap$Population, haplo = h))
(Y <- haploFreq(ali, fac = hap$Year, haplo = h))
(H <- haploFreq(ali, fac = hap$haplotype, haplo = h))


d <- dist.dna(h, "N")
nt <- haploNet(h)

plot(nt)
(sz <- summary(h))
(nt.labs <- attr(nt, "labels"))
sz <- sz[nt.labs]
sz_prop <- sz/sum(sz)

#plot by pop
P <- P[nt.labs, ]
Y <- Y[nt.labs, ]
H <- H[nt.labs, ]


setHaploNetOptions(link.width=1,link.width.alt=0.0,mutations.sequence.length=0.01,mutations.sequence.width=0.5)# this does not wrk properly


pdf("erg24_hap_networks.pdf")
par(mfrow=c(1,2))
gbg <- c("#984EA3","#377EB8","#EA9999","#E41A1C","#E5B110")

plot(nt,xy=xy,size=sz_prop*10,labels=FALSE,pie = P,scale.ratio=1.5,bg=gbg,threshold = c(0,1))

#FLF FVF YLF YLH YVF YVH
gbg <- c("darkgreen", "gold","darkorange3", "blue", "darkred", "darkblue" )


plot(nt,xy=xy,size=sz_prop*10,labels=FALSE,pie = H,scale.ratio=1.5,bg=gbg,threshold = c(0,1))

dev.off()
plot(nt)


xy<- replot() # save coordinates and change manually, the order is the same of h

xy$x[2] <- -5
xy$x[2]
xy$y[2] <- 2.5
#save(xy,file="hap_net_xy.Rdata")

