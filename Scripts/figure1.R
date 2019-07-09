# Figure 1 script.

# Load ribohits data.
ribohits.data <- read.table("Data/ribohits_data_7-8-19.tsv", header = T, sep = "\t", stringsAsFactors = F)

# Load packages.
library(tidyverse)

# Gather length data by frame.
names(ribohits.data)[10:12] <- c("In-frame", "-1 shifted", "-2 shifted")
ribohits.frame <- ribohits.data %>% gather("In-frame", "-1 shifted", "-2 shifted", key = "Frame", value = "Length")

# Group by fragile proteins.
fragile.data <- group_by(ribohits.frame, Fragile, Frame) %>%
  summarise(UTR.mean = mean(log(Length.3UTR), na.rm = T),
            UTR.sd = sd(log(Length.3UTR), na.rm = T),
            UTR.N = length(log(Length.3UTR)),
            UTR.se = UTR.sd / sqrt(UTR.N),
            Length.mean = mean(log(Length), na.rm = T),
            Length.sd = sd(log(Length), na.rm = T),
            Length.N = length(log(Length)),
            Length.se = Length.sd / sqrt(Length.N),
            ISD.Last10.mean = mean(sqrt(ISD.Last10.iupred2), na.rm = T),
            ISD.Last10.sd = sd(sqrt(ISD.Last10.iupred2), na.rm = T),
            ISD.Last10.N = length(sqrt(ISD.Last10.iupred2)),
            ISD.Last10.se = ISD.Last10.sd / sqrt(ISD.Last10.N))
fragile.data$Frame <- factor(fragile.data$Frame)

# Plotting figures.
png(filename = "Scripts/Figures/Fragile_UTRLength_7-8-19.png", width = 350, height = 400)
ggplot(data = fragile.data,
       aes(x = factor(Fragile),
           y = exp(UTR.mean))) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = exp(UTR.mean - UTR.se), ymax = exp(UTR.mean + UTR.se)), width = 0.4, size = 2) +
  scale_x_discrete(labels = c("Not fragile", "Fragile")) +
  xlab("") +
  ylab("3' UTR Length (nts)") +
  scale_y_continuous(limits = c(70, 150)) +
  theme_bw(base_size = 28)
dev.off()

png(filename = "Scripts/Figures/Fragile_ExtL_7-8-19.png", width = 400, height = 400)
ggplot(data = fragile.data,
       aes(x = factor(Fragile),
           y = exp(Length.mean),
           color = Frame)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = exp(Length.mean - Length.se), ymax = exp(Length.mean + Length.se)), width = 0.4, size = 2) +
  scale_x_discrete(labels = c("Not fragile", "Fragile")) +
  xlab("") +
  ylab("Extension length (AAs)") +
  theme_bw(base_size = 28) +
  scale_color_manual(breaks = c("In-frame", "-1 shifted", "-2 shifted"),
                     values = c("#E69F00", "#D55E00", "#009E73")) +
  theme(legend.position = c(0.74, 0.7675),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black"),
        legend.spacing.y = unit(0.25, 'cm'), legend.key.size = unit(1.0, "cm"))
dev.off()

png(filename = "Scripts/Figures/Fragile_ISDLast10_7-8-19.png", width = 350, height = 400)
ggplot(data = fragile.data,
       aes(x = factor(Fragile),
           y = ISD.Last10.mean ^ 2)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = (ISD.Last10.mean - ISD.Last10.se)^2, ymax = (ISD.Last10.mean + ISD.Last10.se)^2),
                width = 0.4, size = 2) +
  scale_x_discrete(labels = c("Not fragile", "Fragile")) +
  scale_y_continuous(limits = c(0.22, 0.29)) +
  xlab("") +
  ylab("C-terminal ISD") +
  theme_bw(base_size = 28)
dev.off()
