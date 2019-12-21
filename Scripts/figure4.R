# Figure 4 code.

# Load ribohits data.
ribohits.data <- read.table("Data/ribohits_data_11-13-19.tab", header = T, sep = "\t", stringsAsFactors = F)

# Load packages.
library(tidyverse)
library(Hmisc)

# Turning ribohits into a factor.
ribohits.data$Ribohits.factor <- ifelse(ribohits.data$Ribohits.Binary == 1, "Yes", "No")
ribohits.data$Ribohits.factor <- factor(ribohits.data$Ribohits.factor, levels = c("Yes", "No"))

# Transforming ISD in frame.
ribohits.data$ISD.0.sqrt <- sqrt(ribohits.data$ISD.0.iupred2)

# Gather into ISD groups.
ribo.isd <-
  ribohits.data[ribohits.data$Fragile == 0,] %>%
  gather("ISD.0.sqrt", "ISD.grafted.sqrt", key = "Grafting", value = "ISD")

isd.data.weighted <- dplyr::group_by(ribo.isd, Grafting, Ribohits.factor) %>%
  dplyr::summarise(ISD.mean = wtd.mean(ISD, weights = Length.0, na.rm = T),
            ISD.sd = sqrt(wtd.var(ISD, weights = Length.0, na.rm = T)),
            ISD.N = length(ISD[!is.na(ISD)]),
            ISD.se = ISD.sd / sqrt(ISD.N))
isd.data.weighted$Grafting <- c("Original", "Original", "Grafted", "Grafted")
isd.data.weighted$Grafting <- factor(isd.data.weighted$Grafting, levels = c("Original", "Grafted"))
names(isd.data.weighted)[2] <- "Ribohits"

# Making figure 3 part A.
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

fig.4a.orig <- 
  ggplot(data = ribohits.data[ribohits.data$Fragile == 0,],
         aes(x = log10(Ribohits.summed + 0.5),
             y = ISD.0.sqrt,
             weight = Length.0 - 1)
  ) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("ISD") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = sqrt(c(0, 0.1, 0.25, 0.5, 1)),
                     labels = c("0.00", "0.10", "0.25", "0.50", "1.00")) +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28)

#ggExtra::ggMarginal(fig.4a.orig, type = "histogram")
png(filename = "Scripts/Figures/ISD_0_weighted_loess_lm_12-8-19.png", width = 6, height = 6, units = "in", res = 350)
fig.4a.orig
dev.off()

# Part A grafted.
png(filename = "Scripts/Figures/ISD_g_weighted_loess_lm_12-8-19.png", width = 6, height = 6, units = "in", res = 350)
fig.4a.grafted <-
  ggplot(data = ribohits.data[ribohits.data$Fragile == 0,],
       aes(x = log10(Ribohits.summed + 0.5),
           y = ISD.grafted.sqrt,
           weight = Length.0)
) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("ISD") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = sqrt(c(0, 0.1, 0.25, 0.5, 1)),
                     labels = c("0.00", "0.10", "0.25", "0.50", "1.00")) +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28)
#ggExtra::ggMarginal(fig.4a.grafted, type = "histogram")
fig.4a.grafted
dev.off()

# Making figure 3 part B.
names(ribohits.data)[19] <- "Ribohits"
isd.0.lm <- lm(data = ribohits.data[ribohits.data$Fragile == 0,],
               formula = ISD.0.sqrt ~ log10(Ribohits.summed + 0.5) + log10(Length.0),
               weights = (Length.0 - 1))
isd.0.summary <- summary(isd.0.lm)
AIC(isd.0.lm)
isd.0.lrt <- drop1(isd.0.lm, test = "Chisq")
isd.0.lrt$`Pr(>Chi)`

isd.g.lm <- lm(data = ribohits.data[ribohits.data$Fragile == 0,],
               formula = ISD.grafted.sqrt ~ log10(Ribohits.summed + 0.5) + log10(Length.0),
               weights = Length.0 - 1)
isd.g.summary <- summary(isd.g.lm)
AIC(isd.g.lm)
isd.g.lrt <- drop1(isd.g.lm, test = "Chisq")
isd.g.lrt$`Pr(>Chi)`

med.ribo <- log10(median(ribohits.data[ribohits.data$Fragile == 0,]$Ribohits.summed + 0.5, na.rm = T))
no.ribo <- log10(0.5)

png(filename = "Scripts/Figures/ISD0_Length0_ribo_12-8-19.png", width = 6, height = 6, units = "in", res = 350)
ggplot(data = ribohits.data[ribohits.data$Fragile == 0,],
       aes(x = log10(Length.0),
           y = sqrt(ISD.0.iupred2),
           color = Ribohits)) +
  stat_function(fun = function(x)isd.0.lm$coefficients[1]+isd.0.lm$coefficients[3]*x+no.ribo*isd.0.lm$coefficients[2],
                geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)isd.0.lm$coefficients[1]+isd.0.lm$coefficients[3]*x+med.ribo*isd.0.lm$coefficients[2],
                geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("ISD") +
  scale_x_continuous(breaks = log10(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  scale_y_continuous(breaks = sqrt(c(0.04, 0.16, 0.36)),
                     labels = c(0.04, 0.16, 0.36),
                     limits = c(0.12, 0.62)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), labels = c("Median", "None"), values = cbPalette[c(2,4)]) +
  theme(legend.position = c(0.79, 0.8485),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black"))
dev.off()

png(filename = "Scripts/Figures/ISDg_Length0_ribo_12-8-19.png", width = 6, height = 6, units = "in", res = 350)
ggplot(data = ribohits.data[ribohits.data$Fragile == 0,],
       aes(x = log10(Length.0),
           y = ISD.grafted.sqrt,
           color = Ribohits)) +
  stat_function(fun = function(x)isd.g.lm$coefficients[1]+isd.g.lm$coefficients[3]*x+no.ribo*isd.0.lm$coefficients[2],
                geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)isd.g.lm$coefficients[1]+isd.g.lm$coefficients[3]*x+med.ribo*isd.g.lm$coefficients[2],
                geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("ISD") +
  scale_x_continuous(breaks = log10(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  scale_y_continuous(breaks = sqrt(c(0.04, 0.16, 0.36)),
                     labels = c(0.04, 0.16, 0.36),
                     limits = c(0.12, 0.62)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), labels = c("Median", "None"), values = cbPalette[c(2,4)]) +
  theme(legend.position = "none")
dev.off()
