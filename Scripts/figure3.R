# Figure 3 script.

# Load ribohits data.
ribohits.data <- read.table("Data/ribohits_data_11-13-19.tab", header = T, sep = "\t", stringsAsFactors = F)

# Load packages.
library(tidyverse)

# Gather length data by frame.
ribohits.data$Fragile.factor <- ifelse(ribohits.data$Fragile == 0, "Backed up", "Lacks backup")
ribohits.data$Fragile.factor <- factor(ribohits.data$Fragile.factor, levels = c("Backed up", "Lacks backup"))
names(ribohits.data)[c(10, 21:22)] <- c("In-frame", "+1 shifted", "+2 shifted")
ribohits.frame <- ribohits.data %>% gather("In-frame", "+1 shifted", "+2 shifted", key = "Frame", value = "Length")

# Group by fragile proteins.
fragile.data <- group_by(ribohits.frame, Fragile, Frame) %>%
  summarise(UTR.mean = mean(log(Length.3UTR), na.rm = T),
            UTR.sd = sd(log(Length.3UTR), na.rm = T),
            UTR.N = length(log(Length.3UTR)),
            UTR.se = UTR.sd / sqrt(UTR.N),
            Length.mean = mean(log(Length), na.rm = T),
            Length.sd = sd(log(Length), na.rm = T),
            Length.N = sum(!is.na(Length)),
            Length.se = Length.sd / sqrt(Length.N),
            ISD.Last10.mean = mean(sqrt(ISD.Last10.iupred2), na.rm = T),
            ISD.Last10.sd = sd(sqrt(ISD.Last10.iupred2), na.rm = T),
            ISD.Last10.N = length(sqrt(ISD.Last10.iupred2)),
            ISD.Last10.se = ISD.Last10.sd / sqrt(ISD.Last10.N),
            tAI.mean = mean(tai.4.log, na.rm = T),
            tAI.sd = sd(tai.4.log, na.rm = T),
            tAI.N = sum(!is.na(tai.4.log)),
            tAI.se = tAI.sd / sqrt(tAI.N))
fragile.data$Frame <- factor(fragile.data$Frame, levels = c("In-frame", "+1 shifted", "+2 shifted"))

# Plotting figures. Means and standard errors for each figure first, then with a boxplot.
# Part A.
png(filename = "Scripts/Figures/Fragile_UTRLength_12-18-19.png", width = 6, height = 6, units = "in", res = 350)
ggplot(data = fragile.data,
       aes(x = factor(Fragile),
           y = exp(UTR.mean))) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = exp(UTR.mean - UTR.se), ymax = exp(UTR.mean + UTR.se)), width = 0.4, size = 2) +
  scale_x_discrete(labels = c("Backed up", "Lacks backup")) +
  xlab("") +
  ylab("3' UTR Length (nts)") +
  scale_y_continuous(limits = c(70, 150)) +
  theme_bw(base_size = 28)
dev.off()

# Part A with box plot.
boxplot.quantiles <- function(x){
  qntls <- quantile(x, probs = c(0.09, 0.25, 0.5, 0.75, 0.91))
  names(qntls) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(qntls)
}
png(filename = "Scripts/Figures/Fragile_UTRLength_boxplot_10-15-19.png", width = 400, height = 400)
ggplot(data = ribohits.data,
       aes(x = factor(Fragile),
           y = log(Length.3UTR))) +
  stat_summary(fun.data = boxplot.quantiles, geom = "boxplot") +
  scale_x_discrete(labels = c("Backed up", "Lacks backup")) +
  xlab("") +
  ylab("3' UTR Length (nts)") +
  scale_y_continuous(breaks = log(c(50, 100, 200)),
                     labels = c(50, 100, 200)) +
  theme_bw(base_size = 28)
dev.off()

# Part B
png(filename = "Scripts/Figures/Fragile_ExtL_+1+2_12-18-19.png", width = 6, height = 6, units = "in", res = 350)
ggplot(data = fragile.data,
       aes(x = factor(Fragile),
           y = exp(Length.mean),
           color = Frame)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = exp(Length.mean - Length.se), ymax = exp(Length.mean + Length.se)), width = 0.4, size = 2) +
  scale_x_discrete(labels = c("Backed up", "Lacks backup")) +
  xlab("") +
  ylab("Extension length (AAs)") +
  theme_bw(base_size = 28) +
  scale_color_manual(breaks = c("In-frame", "+1 shifted", "+2 shifted"),
                     values = c("#E69F00", "#D55E00", "#009E73")) +
  theme(legend.position = c(0.7575, 0.7895),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black"),
        legend.spacing.y = unit(0.25, 'cm'), legend.key.size = unit(1.0, "cm"))
dev.off()

# Part B with box plot.
png(filename = "Scripts/Figures/Fragile_ExtL_boxplot_+1+2_10-27-19.png", width = 600, height = 400)
ggplot(data = ribohits.frame[!is.na(ribohits.frame$Length),],
       aes(x = factor(Fragile),
           y = log(Length),
           fill = Frame)) +
  stat_summary(fun.data = boxplot.quantiles, geom = "boxplot", position = "dodge2") +
  scale_x_discrete(labels = c("Backed up", "Lacks backup")) +
  xlab("") +
  ylab("Extension length (AAs)") +
  scale_y_continuous(breaks = log(c(1, 2, 5, 10, 25)),
                     labels = c(1, 2, 5, 10, 25)) +
  theme_bw(base_size = 28) +
  scale_fill_manual(breaks = c("In-frame", "+1 shifted", "+2 shifted"),
                     values = c("#E69F00", "#D55E00", "#009E73")) +
  theme(legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black"),
        legend.spacing.y = unit(0.25, 'cm'), legend.key.size = unit(1.0, "cm"))
dev.off()

# Part C
png(filename = "Scripts/Figures/Fragile_ISDLast10_12-18-19.png", width = 6, height = 6, units = "in", res = 350)
ggplot(data = fragile.data,
       aes(x = factor(Fragile),
           y = ISD.Last10.mean ^ 2)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = (ISD.Last10.mean - ISD.Last10.se)^2, ymax = (ISD.Last10.mean + ISD.Last10.se)^2),
                width = 0.4, size = 2) +
  scale_x_discrete(labels = c("Backed up", "Lacks backup")) +
  scale_y_continuous(limits = c(0.22, 0.29)) +
  xlab("") +
  ylab("C-terminal ISD") +
  theme_bw(base_size = 28)
dev.off()

# Part C with boxplot
png(filename = "Scripts/Figures/Fragile_ISDLast10_boxplot_10-15-19.png", width = 400, height = 400)
ggplot(data = ribohits.data,
       aes(x = factor(Fragile),
           y = sqrt(ISD.Last10.iupred2))) +
  stat_summary(fun.data = boxplot.quantiles, geom = "boxplot") +
  scale_x_discrete(labels = c("Backed up", "Lacks backup")) +
  scale_y_continuous(breaks = sqrt(c(0.05, 0.15, 0.25, 0.4, 0.6)),
                     labels = c("0.05", "0.15", "0.25", "0.40", "0.60")) +
  xlab("") +
  ylab("C-terminal ISD") +
  theme_bw(base_size = 28)
dev.off()

# Part D
png(filename = "Scripts/Figures/Fragile_tAI4log_12-18-19.png", width = 6, height = 6, units = "in", res = 350)
ggplot(fragile.data[fragile.data$Frame == "In-frame",],
       aes(
         x = factor(Fragile),
         y = exp(tAI.mean)
       )) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = exp(tAI.mean - tAI.se), ymax = exp(tAI.mean + tAI.se)),
                width = 0.4, size = 2) +
  ylab("tAI") +
  xlab("") +
  scale_y_continuous(limits = c(0.278, 0.312)) +
  scale_x_discrete(labels = c("Backed up", "Lacks backup")) +
  theme_bw(base_size = 28)
dev.off()

# Part D with boxplot.
png(filename = "Scripts/Figures/Fragile_tAI4log_boxplot_10-15-19.png", height = 400, width = 425)
ggplot(ribohits.data,
       aes(
         x = factor(Fragile),
         y = tai.4.log
       )) +
  stat_summary(fun.data = boxplot.quantiles, geom = "boxplot") +
  ylab("tAI") +
  xlab("") +
  scale_y_continuous(breaks = log(c(0.2, 0.3, 0.4)),
                     labels = c(0.2, 0.3, 0.4)) +
  scale_x_discrete(labels = c("Backed up", "Lacks backup")) +
  theme_bw(base_size = 28)
dev.off()
