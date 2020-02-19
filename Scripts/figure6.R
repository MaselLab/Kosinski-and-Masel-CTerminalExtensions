# Figure 6 code.

# Load ribohits data.
ribohits.data <- read.table("ribohits_data_11-13-19.tab", header = T, sep = "\t", stringsAsFactors = F)

# Load packages.
library(tidyverse)

# Removing fragile proteins.
ribohits.nofragile <- ribohits.data[ribohits.data$Fragile == 0,]

# Building the disorder propensity and RSA models.
disprop.lm <- lm(
  data = ribohits.nofragile,
  formula = Disorder.prop.0 ~ log10(Ribohits.summed + 0.5) + log10(Length.0),
  weights = Length.0 - 1
)
disprop.summary <- summary(disprop.lm)
disprop.int <- disprop.summary$coefficients[1,1]
disprop.b.ribo <- disprop.summary$coefficients[2,1]
disprop.b.length <- disprop.summary$coefficients[3,1]
drop1(disprop.lm, test = "Chisq")

hydrophilicity.lm <- lm(
  data = ribohits.nofragile,
  formula = RSA.hydrophilicity.0 ~ log10(Ribohits.summed + 0.5) + log10(Length.0),
  weights = Length.0 - 1
)
rsa.summary <- summary(hydrophilicity.lm)
rsa.int <- rsa.summary$coefficients[1,1]
rsa.b.ribo <- rsa.summary$coefficients[2,1]
rsa.b.length <- rsa.summary$coefficients[3,1]
drop1(hydrophilicity.lm, test = "Chisq")

ribo.median <- log10(median(ribohits.nofragile$Ribohits.summed + 0.5, na.rm = T))
ribo.none <- log10(0.5)

# Weighted R2.
disprop.nolength.lm <- lm(
  data = ribohits.nofragile,
  formula = Disorder.prop.0 ~ log10(Ribohits.summed + 0.5),
  weights = Length.0 - 1
)
disprop.nolength.summary <- summary(disprop.nolength.lm)
disprop.nolength.summary

disprop.nolength.noweights.lm <- lm(
  data = ribohits.nofragile,
  formula = Disorder.prop.0 ~ log10(Ribohits.summed + 0.5)
)
disprop.nolength.noweights.summary <- summary(disprop.nolength.noweights.lm)
disprop.nolength.noweights.summary

disprop.noweights.lm <- lm(
  data = ribohits.nofragile,
  formula = Disorder.prop.0 ~ log10(Ribohits.summed + 0.5) + log10(Length.0)
)
disprop.noweights.summary <- summary(disprop.noweights.lm)
disprop.noweights.summary

hydrophilicity.nolength.lm <- lm(
  data = ribohits.nofragile,
  formula = RSA.hydrophilicity.0 ~ log10(Ribohits.summed + 0.5),
  weights = Length.0 - 1
)
rsa.nolength.summary <- summary(hydrophilicity.nolength.lm)
rsa.nolength.summary

hydrophilicity.nolength.noweights.lm <- lm(
  data = ribohits.nofragile,
  formula = RSA.hydrophilicity.0 ~ log10(Ribohits.summed + 0.5)
)
rsa.nolength.noweights.summary <- summary(hydrophilicity.nolength.noweights.lm)
rsa.nolength.noweights.summary

hydrophilicity.noweights.lm <- lm(
  data = ribohits.nofragile,
  formula = RSA.hydrophilicity.0 ~ log10(Ribohits.summed + 0.5) + log10(Length.0)
)
rsa.noweights.summary <- summary(hydrophilicity.noweights.lm)
rsa.noweights.summary

# Part A.
fig.6a.dprop <-
  ggplot(data = ribohits.nofragile,
         aes(x = log10(Ribohits.summed + 0.5),
             y = Disorder.prop.0,
             weight = Length.0 - 1,
             size = Length.0 - 1)
  ) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", color = "#33CCFF", lwd = 1.1) +
  geom_smooth(method = "lm", color = "#CC0000", lwd = 1.1) +
  xlab("Ribohits") +
  ylab("Disorder propensity") +
  #scale_y_continuous(breaks = sqrt(c(0, 0.1, 0.25, 0.5)),
  #                   labels = c("0.00", "0.10", "0.25", "0.50")) +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
png(filename = "Disprop_weighted_loess_lm_2-18-20.png", width = 6, height = 6, units = "in", res = 350)
#ggExtra::ggMarginal(fig.6a.dprop, type = "histogram")
fig.6a.dprop
dev.off()

fig.6a.rsa <-
  ggplot(data = ribohits.nofragile,
         aes(x = log10(Ribohits.summed + 0.5),
             y = RSA.hydrophilicity.0,
             weight = Length.0 - 1,
             size = Length.0 - 1)
  ) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", color = "#33CCFF", lwd = 1.1) +
  geom_smooth(method = "lm", color = "#CC0000", lwd = 1.1) +
  xlab("Ribohits") +
  ylab("RSA") +
  #scale_y_continuous(breaks = sqrt(c(0, 0.1, 0.25, 0.5)),
  #                   labels = c("0.00", "0.10", "0.25", "0.50")) +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
png(filename = "RSA_weighted_loess_lm_2-18-20.png", width = 6, height = 6, units = "in", res = 350)
#ggExtra::ggMarginal(fig.6a.rsa, type = "histogram")
fig.6a.rsa
dev.off()

# Unweighted versions of part A.
fig.6a.unweighted.dprop <-
  ggplot(data = ribohits.nofragile,
         aes(x = log10(Ribohits.summed + 0.5),
             y = Disorder.prop.0)
  ) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("Disorder propensity") +
  #scale_y_continuous(breaks = sqrt(c(0, 0.1, 0.25, 0.5)),
  #                   labels = c("0.00", "0.10", "0.25", "0.50")) +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28)
png(filename = "Disprop_unweighted_loess_lm_2-18-20.png", width = 6, height = 6, units = "in", res = 350)
#ggExtra::ggMarginal(fig.6a.dprop, type = "histogram")
fig.6a.unweighted.dprop
dev.off()

fig.6a.unweighted.rsa <-
  ggplot(data = ribohits.nofragile,
         aes(x = log10(Ribohits.summed + 0.5),
             y = RSA.hydrophilicity.0)
  ) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("RSA") +
  #scale_y_continuous(breaks = sqrt(c(0, 0.1, 0.25, 0.5)),
  #                   labels = c("0.00", "0.10", "0.25", "0.50")) +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28)
png(filename = "RSA_unweighted_loess_lm_2-18-20.png", width = 6, height = 6, units = "in", res = 350)
#ggExtra::ggMarginal(fig.6a.rsa, type = "histogram")
fig.6a.unweighted.rsa
dev.off()

# Part B
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

png(filename = "DisProp_L0_ribo_2-19-20.png", width = 6, height = 6, units = "in", res = 350)
ggplot(data = ribohits.data[ribohits.data$Fragile == 0,],
       aes(x = log10(Length.0),
           y = Disorder.prop.0,
           color = Ribohits)) +
  stat_function(fun = function(x)disprop.int+disprop.b.length*x+disprop.b.ribo*ribo.none,
                geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)disprop.int+disprop.b.length*x+disprop.b.ribo*ribo.median,
                geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("Disorder propensity") +
  scale_x_continuous(breaks = log10(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  #scale_y_continuous(breaks = sqrt(c(0.04, 0.16, 0.36)),
  #                   labels = c(0.04, 0.16, 0.36),
  #                   limits = c(0.125, 0.725)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), values = cbPalette[c(2,6)]) +
  theme(legend.position = "none")
dev.off()

png(filename = "Hydrophilic_L0_ribo_2-19-20.png", width = 6, height = 6, units = "in", res = 350)
ggplot(data = ribohits.data[ribohits.data$Fragile == 0,],
       aes(x = log10(Length.0),
           y = RSA.hydrophilicity.0,
           color = Ribohits)) +
  stat_function(fun = function(x)rsa.int+rsa.b.length*x+rsa.b.ribo*ribo.none,
                geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)rsa.int+rsa.b.length*x+rsa.b.ribo*ribo.median,
                geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("RSA") +
  scale_x_continuous(breaks = log10(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  #scale_y_continuous(breaks = sqrt(c(0.04, 0.16, 0.36)),
  #                   labels = c(0.04, 0.16, 0.36),
  #                   limits = c(0.125, 0.725)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), labels = c("Median", "Zero"), values = cbPalette[c(2,6)]) +
  theme(legend.position = c(0.785, 0.835),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black"))
dev.off()
