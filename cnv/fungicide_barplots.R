# Libraries
library(ggplot2)
library(ggpubr)
library(dplyr)

setwd("~/projects/nikos/fungicide/resistance/")

# cytb

# read metadata
cytB_metadata <- read.csv("2_output/cytB_metadata_new_samples.csv")
cytB_metadata_all_samples <- read.csv("2_output/cytB_metadata.csv")

# exclude Bgh isolates
cytB_metadata_all_samples <- cytB_metadata_all_samples %>% filter(!fs_level_4 == "")

# population
p <- ggplot(cytB_metadata, aes(x=fs_level_4, fill=cytb_G143A)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# country
c <- ggplot(cytB_metadata, aes(x=Country, fill=cytb_G143A)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

countries <- c("United Kingdom", "France", "Switzerland")

cytB_metadata_F_CH_UK <- cytB_metadata_all_samples %>%
  filter(Country %in% countries)

# year of collection
y <- ggplot(cytB_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=cytb_G143A)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# cytb plot
pdf("2_output/cytb_barplots.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, p, ncol = 2), c,
          nrow = 2)
dev.off()

###

# erg24_D137E

# read metadata
erg24_metadata <- read.csv("2_output/erg24_metadata_new_samples.csv")
erg24_metadata_all_samples <- read.csv("2_output/erg24_metadata.csv")

# population
p <- ggplot(erg24_metadata, aes(x=fs_level_4, fill=erg24_D137E)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# country
c <- ggplot(erg24_metadata, aes(x=Country, fill=erg24_D137E)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

countries <- c("United Kingdom", "France", "Switzerland")

erg24_metadata_F_CH_UK <- erg24_metadata_all_samples %>%
  filter(Country %in% countries)

# year of collection
y <- ggplot(erg24_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=erg24_D137E)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# erg24_D137E plot
pdf("2_output/erg24_D137E_barplots_nolegend.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, p, ncol = 2), c,
          nrow = 2)
dev.off()

# erg24_Y165F

# population
p <- ggplot(erg24_metadata, aes(x=fs_level_4, fill=erg24_Y165F)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# country
c <- ggplot(erg24_metadata, aes(x=Country, fill=erg24_Y165F)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

countries <- c("United Kingdom", "France", "Switzerland")

erg24_metadata_F_CH_UK <- erg24_metadata_all_samples %>%
  filter(Country %in% countries)

# year of collection
y <- ggplot(erg24_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=erg24_Y165F)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# erg24_Y165F plot
pdf("2_output/erg24_Y165F_barplots_nolegend.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, p, ncol = 2), c,
          nrow = 2)
dev.off()

# erg24_F289H

# population
p <- ggplot(erg24_metadata, aes(x=fs_level_4, fill=erg24_F289H)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# country
c <- ggplot(erg24_metadata, aes(x=Country, fill=erg24_F289H)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

countries <- c("United Kingdom", "France", "Switzerland")

erg24_metadata_F_CH_UK <- erg24_metadata_all_samples %>%
  filter(Country %in% countries)

# year of collection
y <- ggplot(erg24_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=erg24_F289H)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# erg24_F289H plot
pdf("2_output/erg24_F289H_barplots_nolegend.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, p, ncol = 2), c,
          nrow = 2)
dev.off()

# erg24_D291N

# population
p <- ggplot(erg24_metadata, aes(x=fs_level_4, fill=erg24_D291N)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# country
c <- ggplot(erg24_metadata, aes(x=Country, fill=erg24_D291N)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

countries <- c("United Kingdom", "France", "Switzerland")

erg24_metadata_F_CH_UK <- erg24_metadata_all_samples %>%
  filter(Country %in% countries)

# year of collection
y <- ggplot(erg24_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=erg24_D291N)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# erg24_D291N plot
pdf("2_output/erg24_D291N_barplots_nolegend.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, p, ncol = 2), c,
          nrow = 2)
dev.off()

# erg24_V295L

# population
p <- ggplot(erg24_metadata, aes(x=fs_level_4, fill=erg24_V295L)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577", "black")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# country
c <- ggplot(erg24_metadata, aes(x=Country, fill=erg24_V295L)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577", "black")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

countries <- c("United Kingdom", "France", "Switzerland")

erg24_metadata_F_CH_UK <- erg24_metadata_all_samples %>%
  filter(Country %in% countries)

# year of collection
y <- ggplot(erg24_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=erg24_V295L)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577", "black")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# erg24_V295L plot
pdf("2_output/erg24_V295L_barplots_nolegend.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, p, ncol = 2), c,
          nrow = 2)
dev.off()

###

# cyp51

# read metadata
cyp51_metadata <- read.csv("2_output/cyp51_metadata_new_samples_new_ratio.csv")
cyp51_metadata_all_samples <- read.csv("2_output/cyp51_metadata_new_ratio.csv")
cyp51_metadata_all_samples <- cyp51_metadata_all_samples %>% filter(Sample.Name!="GBR_WC1110")
cyp51_metadata_F_CH_UK <- cyp51_metadata_all_samples %>%
  filter(Country %in% countries)

# S79T population
p <- ggplot(cyp51_metadata, aes(x=fs_level_4, fill=cyp51_S79T)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# S79T country
c <- ggplot(cyp51_metadata, aes(x=Country, fill=cyp51_S79T)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# S79T year of collection
y <- ggplot(cyp51_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=cyp51_S79T)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

cn <- ggplot(cyp51_metadata, aes(x=as.factor(cn), fill=cyp51_S79T)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("Copy number") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.x=element_text(size=20)) +
  theme(legend.position = "none")

# S79T plot
pdf("2_output/cyp51_S79T_barplots_22_8_2024.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, cn, p, ncol = 3), c,
          nrow = 2)
dev.off()

# Y136F population
p <- ggplot(cyp51_metadata, aes(x=fs_level_4, fill=cyp51_Y136F)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# Y136F country
c <- ggplot(cyp51_metadata, aes(x=Country, fill=cyp51_Y136F)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# Y136F year of collection
y <- ggplot(cyp51_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=cyp51_Y136F)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

cn <- ggplot(cyp51_metadata, aes(x=as.factor(cn), fill=cyp51_Y136F)) +
  geom_bar() +
  scale_x_discrete("Copy number") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.x=element_text(size=20)) +
  theme(legend.position = "none")

# Y136F plot
pdf("2_output/cyp51_Y136F_barplots_22_8_2024.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, cn, p, ncol = 3), c,
          nrow = 2)
dev.off()

# K175N population
p <- ggplot(cyp51_metadata, aes(x=fs_level_4, fill=cyp51_K175N)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# K175N country
c <- ggplot(cyp51_metadata, aes(x=Country, fill=cyp51_K175N)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# K175N year of collection
y <- ggplot(cyp51_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=cyp51_K175N)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

cn <- ggplot(cyp51_metadata, aes(x=as.factor(cn), fill=cyp51_K175N)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("Copy number") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.x=element_text(size=20)) +
  theme(legend.position = "none")

# K175N plot
pdf("2_output/cyp51_K175N_barplots_22_8_2024.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, cn, p, ncol = 3), c,
          nrow = 2)
dev.off()

# S509T population
p <- ggplot(cyp51_metadata, aes(x=fs_level_4, fill=cyp51_S509T)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# S509T country
c <- ggplot(cyp51_metadata, aes(x=Country, fill=cyp51_S509T)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# S509T year of collection
y <- ggplot(cyp51_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=cyp51_S509T)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

cn <- ggplot(cyp51_metadata, aes(x=as.factor(cn), fill=cyp51_S509T)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("Copy number") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.x=element_text(size=20)) +
  theme(legend.position = "none")

# S509T plot
pdf("2_output/cyp51_S509T_barplots_22_8_2024.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, cn, p, ncol = 3), c,
          nrow = 2)
dev.off()

# I3K population
p <- ggplot(cyp51_metadata, aes(x=fs_level_4, fill=cyp51_I3K)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# I3K country
c <- ggplot(cyp51_metadata, aes(x=Country, fill=cyp51_I3K)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# I3K year of collection
y <- ggplot(cyp51_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=cyp51_I3K)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

cn <- ggplot(cyp51_metadata, aes(x=as.factor(cn), fill=cyp51_I3K)) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("Copy number") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#ff8577", "#1a82d2")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.x=element_text(size=20)) +
  theme(legend.position = "none")

# K3I plot
pdf("2_output/cyp51_I3K_barplots_22_8_2024.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, cn, p, ncol = 3), c,
          nrow = 2)
dev.off()

# L236F population
p <- ggplot(cyp51_metadata, aes(x=fs_level_4, fill=cyp51_L236F)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# L236F country
c <- ggplot(cyp51_metadata, aes(x=Country, fill=cyp51_L236F)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# L236F year of collection
y <- ggplot(cyp51_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=cyp51_L236F)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

cn <- ggplot(cyp51_metadata, aes(x=as.factor(cn), fill=cyp51_L236F)) +
  geom_bar() +
  scale_x_discrete("Copy number") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.x=element_text(size=20)) +
  theme(legend.position = "none")

# L236F plot
pdf("2_output/cyp51_L236F_barplots_22_8_2024.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, cn, p, ncol = 3), c,
          nrow = 2)
dev.off()

# T271S population
p <- ggplot(cyp51_metadata, aes(x=fs_level_4, fill=cyp51_T271S)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text = element_text(size = 20),
        legend.title=element_text(size=20)) +
  theme(legend.position = "none")

# T271S country
c <- ggplot(cyp51_metadata, aes(x=Country, fill=cyp51_T271S)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# T271S year of collection
y <- ggplot(cyp51_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=cyp51_T271S)) +
  geom_bar() +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

cn <- ggplot(cyp51_metadata, aes(x=as.factor(cn), fill=cyp51_T271S)) +
  geom_bar() +
  scale_x_discrete("Copy number") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#1a82d2", "#ff8577")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.x=element_text(size=20)) +
  theme(legend.position = "none")

# T271S plot
pdf("2_output/cyp51_T271S_barplots_22_8_2024.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, cn, p, ncol = 3), c,
          nrow = 2)
dev.off()

# CN

library(scales)
show_col(viridis_pal()(7))

# population
p <- ggplot(cyp51_metadata, aes(x=fs_level_4, fill=as.factor(cn))) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("") +
  scale_fill_manual(values = c("#FDE725FF", "#8FD744FF", "#35B779FF", "#21908CFF", "#31688EFF", "#443A83FF", "#440154FF")) +
  theme_classic() +
  theme(axis.text=element_text(size=20)) +
  theme(legend.text=element_text(size = 20),
        legend.title=element_text(size=20)) +
  guides(fill=guide_legend("CN")) +
  theme(legend.position = "none")

# country
c <- ggplot(cyp51_metadata, aes(x=Country, fill=as.factor(cn))) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#FDE725FF", "#8FD744FF", "#35B779FF", "#21908CFF", "#31688EFF", "#443A83FF", "#440154FF")) +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# year of collection
y <- ggplot(cyp51_metadata_F_CH_UK, aes(x=Year.of.Collection.Group, fill=as.factor(cn))) +
  geom_bar(position = position_stack(reverse = TRUE)) +
  scale_x_discrete("") +
  scale_y_continuous("Accession count") +
  scale_fill_manual(values = c("#FDE725FF", "#8FD744FF", "#35B779FF", "#21908CFF", "#31688EFF", "#443A83FF")) +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  theme(legend.position = "none")

# plot
pdf("2_output/cyp51_CN_barplots_22_8_2024.pdf", width = 28, height = 12)
ggarrange(ggarrange(y, p, ncol = 2), c,
          nrow = 2)
dev.off()
