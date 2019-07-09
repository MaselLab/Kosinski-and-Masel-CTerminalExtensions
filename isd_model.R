# ISD model script.

# Load ribohits data.
ribohits.data <- read.table("Data/ribohits_data_7-8-19.tsv", header = T, sep = "\t", stringsAsFactors = F)

# Load packages.
library(tidyverse)

# Excluding fragile proteins.
ribo.nofragile <- ribohits.data[ribohits.data$Fragile == 0,]

# Building ISD models and extracting P-values from likelihood ratio tests.
isd.0.lm <- lm(
  data = ribo.nofragile,
  formula = sqrt(ISD.0.iupred2) ~ Ribohits.Binary +log(Length.0),
  weights = Length.0
)
summary(isd.0.lm)
isd.0.test <- drop1(isd.0.lm, test = "Chisq")
isd.0.test$`Pr(>Chi)`
isd.0.no.intercept <- lm(
  data = ribo.nofragile,
  formula = sqrt(ISD.0.iupred2) ~ Ribohits.Binary + log(Length.0) + 0,
  weights = Length.0
)
isd.0.intercept.lrt <- anova(isd.0.no.intercept, isd.0.lm, test = "LRT")
isd.0.intercept.lrt$`Pr(>Chi)`

isd.g.lm <- lm(
  data = ribo.nofragile,
  formula = ISD.grafted.sqrt ~ Ribohits.Binary + log(Length.0),
  weights = Length.0
)
summary(isd.g.lm)
isd.g.test <- drop1(isd.g.lm, test = "Chisq")
isd.g.test$`Pr(>Chi)`
isd.g.no.intercept <- lm(
  data = ribo.nofragile,
  formula = ISD.grafted.sqrt ~ Ribohits.Binary + log(Length.0) + 0,
  weights = Length.0
)
isd.g.intercept.lrt <- anova(isd.g.no.intercept, isd.g.lm, test = "LRT")
isd.g.intercept.lrt$`Pr(>Chi)`
