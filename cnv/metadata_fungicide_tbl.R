# Create a table with metadata and fungicide resistance genes

library(dplyr)
setwd("~/projects/nikos/fungicide/resistance/")

# read metadata
metadata1 <- read.csv("0_data/general/2022+before2022+2023+ncsu_metadata+fs+admxK9_03052024.csv")

# read cytB gene
cytB <- read.csv("2_output/cytb_amino_variable.csv")
cytB <- head(cytB, -5) # remove the last 5 (not tritici)

# keep important columns
metadata <- metadata1 %>%
  select(Sample.Name, Latitude,
         Longitude, fs_level_4,
         Year.of.Collection, Country)

# create cytB metadata table 
cytB_metadata_tbl <- cbind(metadata[match(cytB$ID, metadata$Sample.Name),], cytB$Position_143)
colnames(cytB_metadata_tbl)[7] <- "cytb_G143A"
cytB_metadata_tbl[cytB_metadata_tbl[1]==("GBR_JIW11"), 5] <- 1985
cytB_metadata_tbl[cytB_metadata_tbl[1]==("JIW2"), 5] <- 1980
cytB_metadata_tbl[cytB_metadata_tbl[1]==("GBR_JIW48"), 5] <- 2001

# delete this one
cytB_metadata_tbl <- cytB_metadata_tbl %>% filter(Sample.Name!="GBR_WC1110")

# new column Sensitive vs Resistant
for (i in 1:nrow(cytB_metadata_tbl)) {
  if (cytB_metadata_tbl$cytb_G143A[i] == "G") {
    cytB_metadata_tbl$QoI_interpretation[i] <- "sensitive"  
  } else {
    cytB_metadata_tbl$QoI_interpretation[i] <- "resistant"
  }
}

# new column Sample Year Group
for (i in 1:nrow(cytB_metadata_tbl)) {
  if (cytB_metadata_tbl$Year.of.Collection[i] > 2021) {
    cytB_metadata_tbl$Year.of.Collection.Group[i] <- "2022-2023"  
  } else if (cytB_metadata_tbl$Year.of.Collection[i] > 2006 & cytB_metadata_tbl$Year.of.Collection[i] < 2022) {
    cytB_metadata_tbl$Year.of.Collection.Group[i] <- "2007-2020"
  } else {
    cytB_metadata_tbl$Year.of.Collection.Group[i] <- "1980-2001"
  }
}

# new table only with new samples
cytB_metadata_tbl_new_samples <- cytB_metadata_tbl %>% filter(Year.of.Collection > 2014)

# save tables
write.csv(cytB_metadata_tbl, "2_output/cytB_metadata.csv", row.names = F)
write.csv(cytB_metadata_tbl_new_samples, "2_output/cytB_metadata_new_samples.csv", row.names = F)

# read cyp51 gene
cyp51 <- read.csv("2_output/cyp51_amino_het_variable.csv")
cyp51 <- cyp51 %>% filter(ID!="GBR_WC1110")

# keep important columns
metadata <- metadata %>%
  select(Sample.Name, Latitude,
         Longitude, fs_level_4,
         Year.of.Collection, Country)

cyp51_cn <- read.csv("0_data/general/cyp51_cnv.csv")
cyp51_cn <- cyp51_cn %>% filter(Isolate!="GBR_WC1110")

# create cyp51 metadata table 

## S79T
cyp51_metadata_tbl <- cbind(metadata[match(cyp51$ID, metadata$Sample.Name),], cyp51$Position_79)
cyp51_metadata_tbl[cyp51_metadata_tbl[1]==("GBR_JIW11"), 5] <- 1985
cyp51_metadata_tbl[cyp51_metadata_tbl[1]==("JIW2"), 5] <- 1980
cyp51_metadata_tbl[cyp51_metadata_tbl[1]==("GBR_JIW48"), 5] <- 2001

# delete this one
cyp51_metadata_tbl <- cyp51_metadata_tbl %>% filter(Sample.Name!="GBR_WC1110")

colnames(cyp51_metadata_tbl)[7] <- "cyp51_S79T"

# new column Sensitive vs Resistant
for (i in 1:nrow(cyp51_metadata_tbl)) {
  if (cyp51_metadata_tbl$cyp51_S79T[i] == "S") {
    cyp51_metadata_tbl$cyp51_S79T_interpretation[i] <- "sensitive"  
  } else {
    cyp51_metadata_tbl$cyp51_S79T_interpretation[i] <- "resistant"
  }
}

## Y136F
cyp51_metadata_tbl <- cbind(cyp51_metadata_tbl[match(cyp51$ID, cyp51_metadata_tbl$Sample.Name),], cyp51$Position_136)
colnames(cyp51_metadata_tbl)[9] <- "cyp51_Y136F"

# new column Sensitive vs Resistant
for (i in 1:nrow(cyp51_metadata_tbl)) {
  if (cyp51_metadata_tbl$cyp51_Y136F[i] == "Y") {
    cyp51_metadata_tbl$cyp51_Y136F_interpretation[i] <- "sensitive"  
  } else {
    cyp51_metadata_tbl$cyp51_Y136F_interpretation[i] <- "resistant"
  }
}

## K175N
cyp51_metadata_tbl <- cbind(cyp51_metadata_tbl[match(cyp51$ID, cyp51_metadata_tbl$Sample.Name),], cyp51$Position_175)
colnames(cyp51_metadata_tbl)[11] <- "cyp51_K175N"

# new column Sensitive vs Resistant
for (i in 1:nrow(cyp51_metadata_tbl)) {
  if (cyp51_metadata_tbl$cyp51_K175N[i] == "K") {
    cyp51_metadata_tbl$cyp51_K175N_interpretation[i] <- "sensitive"  
  } else {
    cyp51_metadata_tbl$cyp51_K175N_interpretation[i] <- "resistant"
  }
}

## S509T
cyp51_metadata_tbl <- cbind(cyp51_metadata_tbl[match(cyp51$ID, cyp51_metadata_tbl$Sample.Name),], cyp51$Position_509)
colnames(cyp51_metadata_tbl)[13] <- "cyp51_S509T"

# new column Sensitive vs Resistant
for (i in 1:nrow(cyp51_metadata_tbl)) {
  if (cyp51_metadata_tbl$cyp51_S509T[i] == "S") {
    cyp51_metadata_tbl$cyp51_S509T_interpretation[i] <- "sensitive"  
  } else {
    cyp51_metadata_tbl$cyp51_S509T_interpretation[i] <- "resistant"
  }
}

## K3I
cyp51_metadata_tbl <- cbind(cyp51_metadata_tbl[match(cyp51$ID, cyp51_metadata_tbl$Sample.Name),], cyp51$Position_3)
colnames(cyp51_metadata_tbl)[15] <- "cyp51_K3I"

# new column Sensitive vs Resistant
for (i in 1:nrow(cyp51_metadata_tbl)) {
  if (cyp51_metadata_tbl$cyp51_K3I[i] == "K") {
    cyp51_metadata_tbl$cyp51_K3I_interpretation[i] <- "sensitive"  
  } else {
    cyp51_metadata_tbl$cyp51_K3I_interpretation[i] <- "resistant"
  }
}

## L236F
cyp51_metadata_tbl <- cbind(cyp51_metadata_tbl[match(cyp51$ID, cyp51_metadata_tbl$Sample.Name),], cyp51$Position_236)
colnames(cyp51_metadata_tbl)[17] <- "cyp51_L236F"

# new column Sensitive vs Resistant
for (i in 1:nrow(cyp51_metadata_tbl)) {
  if (cyp51_metadata_tbl$cyp51_L236F[i] == "L") {
    cyp51_metadata_tbl$cyp51_L236F_interpretation[i] <- "sensitive"  
  } else {
    cyp51_metadata_tbl$cyp51_L236F_interpretation[i] <- "resistant"
  }
}

## T271S
cyp51_metadata_tbl <- cbind(cyp51_metadata_tbl[match(cyp51$ID, cyp51_metadata_tbl$Sample.Name),], cyp51$Position_271)
colnames(cyp51_metadata_tbl)[19] <- "cyp51_T271S"

# new column Sensitive vs Resistant
for (i in 1:nrow(cyp51_metadata_tbl)) {
  if (cyp51_metadata_tbl$cyp51_T271S[i] == "T") {
    cyp51_metadata_tbl$cyp51_T271S_interpretation[i] <- "sensitive"  
  } else {
    cyp51_metadata_tbl$cyp51_T271S_interpretation[i] <- "resistant"
  }
}

# column of counting number of resistant mutations (of the four known)
for (i in 1:nrow(cyp51_metadata_tbl)) {
  count <- 0
  if (cyp51_metadata_tbl$cyp51_S79T_interpretation[i] == "resistant") {
    count <- count + 1
  }
  if (cyp51_metadata_tbl$cyp51_Y136F_interpretation[i] == "resistant") {
    count <- count + 1
  }
  if (cyp51_metadata_tbl$cyp51_K175N_interpretation[i] == "resistant") {
    count <- count + 1
  }
  if (cyp51_metadata_tbl$cyp51_S509T_interpretation[i] == "resistant") {
    count <- count + 1
  }
  cyp51_metadata_tbl$count_of_resistant_mutations[i] <- count
}

# add new columns copy number variation ratio and cn
cyp51_metadata_tbl$cn_ratio <- cyp51_cn$ratio
for (i in 1:nrow(cyp51_metadata_tbl)) {
  if (cyp51_metadata_tbl$cn_ratio[i] < 1.5) {
    cyp51_metadata_tbl$cn[i] <- 1
  }
  if (cyp51_metadata_tbl$cn_ratio[i] > 1.5 & cyp51_metadata_tbl$cn_ratio[i] <= 2.5) {
    cyp51_metadata_tbl$cn[i] <- 2
  }
  if (cyp51_metadata_tbl$cn_ratio[i] > 2.5 & cyp51_metadata_tbl$cn_ratio[i] <= 3.5) {
    cyp51_metadata_tbl$cn[i] <- 3
  }
  if (cyp51_metadata_tbl$cn_ratio[i] > 3.5 & cyp51_metadata_tbl$cn_ratio[i] <= 4.5) {
    cyp51_metadata_tbl$cn[i] <- 4
  }
  if (cyp51_metadata_tbl$cn_ratio[i] > 4.5 & cyp51_metadata_tbl$cn_ratio[i] <= 5.5) {
    cyp51_metadata_tbl$cn[i] <- 5
  }
  if (cyp51_metadata_tbl$cn_ratio[i] > 5.5 & cyp51_metadata_tbl$cn_ratio[i] <= 6.5) {
    cyp51_metadata_tbl$cn[i] <- 6
  }
  if (cyp51_metadata_tbl$cn_ratio[i] > 6.6) {
    cyp51_metadata_tbl$cn[i] <- 7
  }
}

# new column Sample Year Group
for (i in 1:nrow(cyp51_metadata_tbl)) {
  if (cyp51_metadata_tbl$Year.of.Collection[i] > 2021) {
    cyp51_metadata_tbl$Year.of.Collection.Group[i] <- "2022-2023"  
  } else if (cyp51_metadata_tbl$Year.of.Collection[i] > 2006 & cyp51_metadata_tbl$Year.of.Collection[i] < 2022) {
    cyp51_metadata_tbl$Year.of.Collection.Group[i] <- "2007-2020"
  } else {
    cyp51_metadata_tbl$Year.of.Collection.Group[i] <- "1980-2001"
  }
}

# new table only with new samples
cyp51_metadata_tbl_new_samples <- cyp51_metadata_tbl %>% filter(Year.of.Collection > 2014)

# save tables
write.csv(cyp51_metadata_tbl, "2_output/cyp51_metadata.csv", row.names = F)
write.csv(cyp51_metadata_tbl_new_samples, "2_output/cyp51_metadata_new_samples.csv", row.names = F)

# add haplotype_cyp51s to the table
cyp51_tbl <- cyp51_metadata_tbl

head(cyp51_tbl)

# new column Sample Year Group
for (i in 1:nrow(cyp51_tbl)) {
  if (cyp51_tbl$Year.of.Collection[i] > 2014) {
    cyp51_tbl$Year.of.Collection.Group_2[i] <- "2015-2023"
  } else {
    cyp51_tbl$Year.of.Collection.Group_2[i] <- "1980-2010"
  }
}


# only keep columns of interest
cyp51_tbl_aa <- cyp51_tbl

# column of counting number of resistant mutations (of the four known)
for (i in 1:nrow(cyp51_tbl_aa)) {
  if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "resistant" &
      cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "resistant" &
      cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "resistant" &
      cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "resistant") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "h"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "sensitive") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "e"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "resistant") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "f"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "resistant") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "i"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "resistant") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "g"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "resistant") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "m"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "sensitive") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "j"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "sensitive") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "c"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "resistant") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "b"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "sensitive") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "d"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "resistant") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "k"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "sensitive") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "l"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "sensitive") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "a"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "resistant" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "sensitive") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "n"
  }
  else if (cyp51_tbl_aa$cyp51_S79T_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_Y136F_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_K175N_interpretation[i] == "sensitive" &
           cyp51_tbl_aa$cyp51_S509T_interpretation[i] == "resistant") {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "o"
  }
  else {
    cyp51_tbl_aa$haplotype_cyp51[i] <- "p"
  }
}

# final metadata table
matched_cytb_tbl <- cytB_metadata_tbl[match(cytB_metadata_tbl$Sample.Name, cyp51_tbl_aa$Sample.Name),]
final_metadata_tbl <- cbind(matched_cytb_tbl%>%select(!Year.of.Collection.Group),
                            cyp51_tbl_aa[,7:length(cyp51_tbl_aa)])
write_csv(final_metadata_tbl, "final_metadata_tbl.csv")
