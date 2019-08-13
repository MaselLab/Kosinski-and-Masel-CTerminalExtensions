# Disorder propensity and RSA models script.

# Load ribohits data.
ribohits.data <- read.table("Data/ribohits_data_tai_8-12-19.tsv", header = T, sep = "\t", stringsAsFactors = F)

# Load packages.
library(tidyverse)

# Excluding fragile proteins.
ribo.nofragile <- ribohits.data[ribohits.data$Fragile == 0,]

# Disorder propensity and mean RSA model building and likelihood-ratio test P-value extraction.
disprop.0.lm <- lm(
  data = ribo.nofragile,
  formula = Disorder.prop.0 ~ Ribohits.Binary + log(Length.0),
  weights = Length.0
)
summary(disprop.0.lm)
drop1(disprop.0.lm, test = "Chisq")
disprop.0.no.intercept <- lm(
  data = ribo.nofragile,
  formula = Disorder.prop.0 ~ Ribohits.Binary + log(Length.0) + 0,
  weights = Length.0
)
disprop.0.no.intercept.lrt <- anova(disprop.0.no.intercept, disprop.0.lm, test = "LRT")
disprop.0.no.intercept.lrt$`Pr(>Chi)`

RSA.mean.0.lm <- lm(
  data = ribo.nofragile,
  formula = RSA.hydrophilicity.0 ~ Ribohits.Binary + log(Length.0),
  weights = Length.0
)
summary(RSA.mean.0.lm)
drop1(RSA.mean.0.lm, test = "Chisq")
RSA.mean.0.no.intercept <- lm(
  data = ribo.nofragile,
  formula = RSA.hydrophilicity.0 ~ Ribohits.Binary + log(Length.0) + 0,
  weights = Length.0
)
RSA.mean.0.no.intercept.lrt <- anova(RSA.mean.0.no.intercept, RSA.mean.0.lm, test = "LRT")
RSA.mean.0.no.intercept.lrt$`Pr(>Chi)`
