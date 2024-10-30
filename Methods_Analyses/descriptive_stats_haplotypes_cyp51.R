# discriptive stats haplotypes cyp51

library(dplyr)
setwd("~/projects/nikos/fungicide/resistance/2_output/")

cyp51_tbl_aa <- read.csv("Supplementary_Data_S1.csv")

export_cyp51_yg1 <- cyp51_tbl_aa %>%
  group_by(Year.of.Collection.Group, haplotype_cyp51) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = Year.of.Collection.Group,
              values_from = n)

write_csv(export_cyp51_yg1, "cyp51_amino_acid_combinations_year_Group.csv")

export_cyp51_yg2 <- cyp51_tbl_aa %>%
  group_by(Year.of.Collection.Group_2, haplotype_cyp51) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = Year.of.Collection.Group_2,
              values_from = n)

write_csv(export_cyp51_yg2, "cyp51_amino_acid_combinations_year_Group2.csv")
