# t-test or Mann-Whitney U-test for CN in cyp51

setwd("~/projects/nikos/fungicide/resistance/2_output/")
library(tidyverse)

# dataset prep
original_tbl <- read.csv("Supplementary_Data_S1.csv")
head(original_tbl)
countries <- c("United Kingdom", "France", "Switzerland")
new_tbl <- original_tbl %>%
  filter(Country %in% countries) %>%
  filter(Year.of.Collection.Group == "1980-2001" | Year.of.Collection.Group == "2022-2023") %>%
  select(Sample.Name, Year.of.Collection.Group, cn)
class(new_tbl$cn)
new_tbl$Year.of.Collection.Group <- as.factor(new_tbl$Year.of.Collection.Group)

# check normality and variance
qqnorm(new_tbl$cn, pch = 1, frame = FALSE)
qqline(new_tbl$cn, col = "steelblue", lwd = 2)

# run Shapiro-Wilk test
shapiro.test(x = new_tbl$cn)

# H0: The variances of the two samples wt~am are not different
var.test(new_tbl$cn~new_tbl$Year.of.Collection.Group)

# Computing t-test,when assumption of equal variance is violated
tst <- t.test(new_tbl$cn~new_tbl$Year.of.Collection.Group, var.equal = F)
tst

