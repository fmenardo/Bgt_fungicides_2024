# fisher's exact test for temporal_dataset

library(tidyverse)

setwd("~/projects/nikos/fungicide/resistance/2_output/")

my_tbl <- read.csv("Table_S2_temporal_dataset_14_8_24.csv", sep = ",")

head(my_tbl)

my_tbl <- my_tbl %>% filter(Year.of.Collection.Group == "1980-2001" | Year.of.Collection.Group == "2022-2023")

# cytb_G143A
f_cytbG143A_1 <- my_tbl %>%
  group_by(cytb_G143A, Year.of.Collection.Group) %>%
  count()

f_cytbG143A <-
  matrix(c(f_cytbG143A_1$n[1], f_cytbG143A_1$n[3], f_cytbG143A_1$n[2], f_cytbG143A_1$n[4]),
         nrow = 2,
         dimnames = list(allele = c("A", "G"),
                         year = c("old", "new")))

fisher.test(f_cytbG143A, alternative = "t")

# erg24_V295L
f_erg24_V295L_1 <- my_tbl %>%
  group_by(erg24_V295L, Year.of.Collection.Group) %>%
  count()

f_erg24_V295L <-
  matrix(c(f_erg24_V295L_1$n[1], f_erg24_V295L_1$n[3], f_erg24_V295L_1$n[2], f_erg24_V295L_1$n[4]),
         nrow = 2,
         dimnames = list(allele = c("L", "V"),
                         year = c("old", "new")))

fisher.test(f_erg24_V295L, alternative = "t")

# erg24_F289H
f_erg24_F289H_1 <- my_tbl %>%
  group_by(erg24_F289H, Year.of.Collection.Group) %>%
  count()

f_erg24_F289H_1[4,] <- f_erg24_F289H_1[3,]
f_erg24_F289H_1[3,3] <- 0
f_erg24_F289H_1[3,2] <- "1980-2001"

f_erg24_F289H <-
  matrix(c(f_erg24_F289H_1$n[3], f_erg24_F289H_1$n[1], f_erg24_F289H_1$n[4], f_erg24_F289H_1$n[2]),
         nrow = 2,
         dimnames = list(allele = c("H", "F"),
                         year = c("old", "new")))

f_erg24_F289H[is.na(f_erg24_F289H)] <- 0

fisher.test(f_erg24_F289H, alternative = "t")

# erg24_Y165F
f_erg24_Y165F_1 <- my_tbl %>%
  group_by(erg24_Y165F, Year.of.Collection.Group) %>%
  count()

f_erg24_Y165F <-
  matrix(c(f_erg24_Y165F_1$n[1], f_erg24_Y165F_1$n[3], f_erg24_Y165F_1$n[2], f_erg24_Y165F_1$n[4]),
         nrow = 2,
         dimnames = list(allele = c("F", "Y"),
                         year = c("old", "new")))

fisher.test(f_erg24_Y165F, alternative = "t")

# cyp51_Y136F
f_cyp51_Y136F_1 <- my_tbl %>%
  group_by(cyp51_Y136F, Year.of.Collection.Group) %>%
  count()

f_cyp51_Y136F <-
  matrix(c(f_cyp51_Y136F_1$n[1], f_cyp51_Y136F_1$n[3], f_cyp51_Y136F_1$n[2], f_cyp51_Y136F_1$n[4]),
         nrow = 2,
         dimnames = list(allele = c("F", "Y"),
                         year = c("old", "new")))

fisher.test(f_cyp51_Y136F, alternative = "t")

# cyp51_S79T
f_cyp51_S79T_1 <- my_tbl %>%
  group_by(cyp51_S79T, Year.of.Collection.Group) %>%
  count()

f_cyp51_S79T <-
  matrix(c(f_cyp51_S79T_1$n[3], f_cyp51_S79T_1$n[1], f_cyp51_S79T_1$n[4], f_cyp51_S79T_1$n[2]),
         nrow = 2,
         dimnames = list(allele = c("T", "S"),
                         year = c("old", "new")))

fisher.test(f_cyp51_S79T, alternative = "t")

# cyp51_K175N
f_cyp51_K175N_1 <- my_tbl %>%
  group_by(cyp51_K175N, Year.of.Collection.Group) %>%
  count()

f_cyp51_K175N <-
  matrix(c(f_cyp51_K175N_1$n[3], f_cyp51_K175N_1$n[1], f_cyp51_K175N_1$n[4], f_cyp51_K175N_1$n[2]),
         nrow = 2,
         dimnames = list(allele = c("N", "K"),
                         year = c("old", "new")))

fisher.test(f_cyp51_K175N, alternative = "t")

# cyp51_S509T
f_cyp51_S509T_1 <- my_tbl %>%
  group_by(cyp51_S509T, Year.of.Collection.Group) %>%
  count()

f_cyp51_S509T_1[4,] <- f_cyp51_S509T_1[3,]
f_cyp51_S509T_1[3,3] <- 0
f_cyp51_S509T_1[3,2] <- "1980-2001"

f_cyp51_S509T <-
  matrix(c(f_cyp51_S509T_1$n[3], f_cyp51_S509T_1$n[1], f_cyp51_S509T_1$n[4], f_cyp51_S509T_1$n[2]),
         nrow = 2,
         dimnames = list(allele = c("T", "S"),
                         year = c("old", "new")))

f_cyp51_S509T[is.na(f_cyp51_S509T)] <- 0

fisher.test(f_cyp51_S509T, alternative = "t")
