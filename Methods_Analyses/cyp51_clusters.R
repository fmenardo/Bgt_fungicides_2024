setwd("projects")
load("cyp51_clusters.R")
library(stringr)
library(dplyr)

# First use this script.
# Then use the python script "fasta_to_variables_table.py".

# cyp51_het_as_X <- read.csv("nikos/fungicide/resistance/0_data/alignments/het_as_X/Bgt_Eur+_cyp51_p.fa")


inds_cl1 <- my_i_clusters$clusters[[1]]
inds_cl2 <- my_i_clusters$clusters[[2]]
inds_cl3 <- my_i_clusters$clusters[[3]]
inds_cl4 <- my_i_clusters$clusters[[4]]
inds_cl5 <- my_i_clusters$clusters[[5]]
inds_cl6 <- my_i_clusters$clusters[[6]]
inds_cl7 <- my_i_clusters$clusters[[7]]
inds_cl8 <- my_i_clusters$clusters[[8]]

for (i in 1:length(my_i_clusters$clusters[[1]])) {
  inds_cl1[i] <- print(strsplit(my_i_clusters$clusters[[1]][i], "[/]")[[1]][1])
}

write.csv(inds_cl1, "nikos/fungicide/resistance/2_output/cyp51_IBD_cluster1_indlist.txt", row.names = F, quote = F)

for (i in 1:length(my_i_clusters$clusters[[2]])) {
  inds_cl2[i] <- print(strsplit(my_i_clusters$clusters[[2]][i], "[/]")[[1]][1])
}

write.csv(inds_cl2, "nikos/fungicide/resistance/2_output/cyp51_IBD_cluster2_indlist.txt", row.names = F, quote = F)

for (i in 1:length(my_i_clusters$clusters[[3]])) {
  inds_cl3[i] <- print(strsplit(my_i_clusters$clusters[[3]][i], "[/]")[[1]][1])
}

write.csv(inds_cl3, "nikos/fungicide/resistance/2_output/cyp51_IBD_cluster3_indlist.txt", row.names = F, quote = F)

for (i in 1:length(my_i_clusters$clusters[[4]])) {
  inds_cl4[i] <- print(strsplit(my_i_clusters$clusters[[4]][i], "[/]")[[1]][1])
}

write.csv(inds_cl4, "nikos/fungicide/resistance/2_output/cyp51_IBD_cluster4_indlist.txt", row.names = F, quote = F)

for (i in 1:length(my_i_clusters$clusters[[5]])) {
  inds_cl5[i] <- print(strsplit(my_i_clusters$clusters[[5]][i], "[/]")[[1]][1])
}

write.csv(inds_cl5, "nikos/fungicide/resistance/2_output/cyp51_IBD_cluster5_indlist.txt", row.names = F, quote = F)

for (i in 1:length(my_i_clusters$clusters[[6]])) {
  inds_cl6[i] <- print(strsplit(my_i_clusters$clusters[[6]][i], "[/]")[[1]][1])
}

write.csv(inds_cl6, "nikos/fungicide/resistance/2_output/cyp51_IBD_cluster6_indlist.txt", row.names = F, quote = F)

for (i in 1:length(my_i_clusters$clusters[[7]])) {
  inds_cl7[i] <- print(strsplit(my_i_clusters$clusters[[7]][i], "[/]")[[1]][1])
}

write.csv(inds_cl7, "nikos/fungicide/resistance/2_output/cyp51_IBD_cluster7_indlist.txt", row.names = F, quote = F)

for (i in 1:length(my_i_clusters$clusters[[8]])) {
  inds_cl8[i] <- print(strsplit(my_i_clusters$clusters[[8]][i], "[/]")[[1]][1])
}

write.csv(inds_cl8, "nikos/fungicide/resistance/2_output/cyp51_IBD_cluster8_indlist.txt", row.names = F, quote = F)

# simple stats

cluster2_vars <- read.csv("nikos/fungicide/resistance/2_output/cyp51_amino_X_IBDcluster2_variable.csv")
sum(cluster2_vars$Position_509 == "T")
sum(cluster2_vars$Position_509 == "S")
sum(cluster2_vars$Position_509 == "X")
length(cluster2_vars$Position_509)



