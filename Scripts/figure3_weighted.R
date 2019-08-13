# Figure 3 code.

# Load ribohits data.
ribohits.data <- read.table("Data/ribohits_data_tai_8-12-19.tsv", header = T, sep = "\t", stringsAsFactors = F)

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
png(filename = "Scripts/Figures/ISD_0_g_weighted_7-25-19.png", width = 400, height = 400)
ggplot(data = isd.data.weighted,
       aes(x = Grafting,
           y = ISD.mean ^ 2,
           color = Ribohits)
) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = (ISD.mean - ISD.se)^2, ymax = (ISD.mean + ISD.se)^2), width = 0.4, size = 2) +
  xlab("") +
  ylab("ISD") +
  #scale_y_continuous(limits = c(0.075, 0.205)) +
  scale_color_manual(values = cbPalette[c(2,4)]) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

# Making figure 3 part B.
names(ribohits.data)[17] <- "Ribohits"
isd.0.lm <- lm(data = ribohits.data[ribohits.data$Fragile == 0,],
               formula = ISD.0.sqrt ~ Ribohits.Binary + log(Length.0),
               weights = Length.0)
isd.g.lm <- lm(data = ribohits.data[ribohits.data$Fragile == 0,],
               formula = ISD.grafted.sqrt ~ Ribohits.Binary + log(Length.0),
               weights = Length.0)

png(filename = "Scripts/Figures/ISD0_Length0_ribo_7-25-19.png", height = 400, width = 400)
ggplot(data = ribohits.data[ribohits.data$Fragile == 0,],
       aes(x = log(Length.0),
           y = sqrt(ISD.0.iupred2),
           color = Ribohits)) +
  stat_function(fun = function(x)0.53108-0.08208*x, geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)0.53108-0.08208*x+0.04996, geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("ISD") +
  scale_x_continuous(breaks = log(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  scale_y_continuous(breaks = sqrt(c(0.04, 0.16, 0.36)),
                     labels = c(0.04, 0.16, 0.36),
                     limits = c(0.125, 0.62)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), values = cbPalette[c(2,4)]) +
  theme(legend.position = c(0.7775, 0.8325),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black"))
dev.off()

png(filename = "Scripts/Figures/ISDg_Length0_ribo_7-25-19.png", height = 400, width = 400)
ggplot(data = ribohits.data[ribohits.data$Fragile == 0,],
       aes(x = log(Length.0),
           y = ISD.grafted.sqrt,
           color = Ribohits)) +
  stat_function(fun = function(x)0.547272-0.080698*x, geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)0.547272-0.080698*x+0.035126, geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("ISD") +
  scale_x_continuous(breaks = log(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  scale_y_continuous(breaks = sqrt(c(0.04, 0.16, 0.36)),
                     labels = c(0.04, 0.16, 0.36),
                     limits = c(0.125, 0.62)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), values = cbPalette[c(2,4)]) +
  theme(legend.position = "none")
dev.off()
