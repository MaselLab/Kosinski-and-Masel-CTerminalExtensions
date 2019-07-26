# Figure 4 code.

# Load ribohits data.
ribohits.data <- read.table("Data/ribohits_data_7-8-19.tsv", header = T, sep = "\t", stringsAsFactors = F)

# Load packages.
library(tidyverse)

# Removing fragile proteins.
ribohits.nofragile <- ribohits.data[ribohits.data$Fragile == 0,]

# Building the disorder propensity and RSA models.
disprop.lm <- lm(
  data = ribohits.nofragile,
  formula = Disorder.prop.0 ~ Ribohits.Binary + log(Length.0),
  weights = Length.0
)
summary(disprop.lm)

hydrophilicity.lm <- lm(
  data = ribohits.nofragile,
  formula = RSA.hydrophilicity.0 ~ Ribohits.Binary + log(Length.0),
  weights = Length.0
)
summary(hydrophilicity.lm)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

png(filename = "Scripts/Figures/DisProp_L0_ribo_7-25-19.png", height = 400, width = 400)
ggplot(data = ribohits.data[ribohits.data$Fragile == 0,],
       aes(x = log(Length.0),
           y = Disorder.prop.0,
           color = Ribohits)) +
  stat_function(fun = function(x)0.354886-0.006060*x, geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)0.354886-0.006060*x+0.011824, geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("Disorder propensity") +
  scale_x_continuous(breaks = log(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  #scale_y_continuous(breaks = sqrt(c(0.04, 0.16, 0.36)),
  #                   labels = c(0.04, 0.16, 0.36),
  #                   limits = c(0.125, 0.725)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), values = cbPalette[c(2,4)]) +
  theme(legend.position = "none")
dev.off()

png(filename = "Scripts/Figures/Hydrophilic_L0_ribo_7-25-19.png", width = 400, height = 400)
ggplot(data = ribohits.data[ribohits.data$Fragile == 0,],
       aes(x = log(Length.0),
           y = RSA.hydrophilicity.0,
           color = Ribohits)) +
  stat_function(fun = function(x)0.2525676-0.0053003*x, geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)0.2525676-0.0053003*x+0.0053970, geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("RSA") +
  scale_x_continuous(breaks = log(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  #scale_y_continuous(breaks = sqrt(c(0.04, 0.16, 0.36)),
  #                   labels = c(0.04, 0.16, 0.36),
  #                   limits = c(0.125, 0.725)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), values = cbPalette[c(2,4)]) +
  theme(legend.position = c(0.7775, 0.8325),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black"))
dev.off()
