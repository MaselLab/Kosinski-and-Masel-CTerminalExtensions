# Supplemental histograms.

# Load ribohits data.
ribohits.data <- read.table("Data/ribohits_data_11-13-19.tab", header = T, sep = "\t", stringsAsFactors = F)

# Load packages.
library(tidyverse)

# Global variables.
today.date <- "1-24-20"

# Histogram for log-transformed 3' UTR hits.
fragile.3utr.filename <- paste("Scripts/Figures/fragile_3utr_", today.date, ".png", sep = "")
png(fragile.3utr.filename, height = 500, width = 500)
ggplot(
  data = ribohits.data,
  aes(
    x = log(Length.3UTR),
    fill = factor(Fragile)
  )
) +
  geom_histogram(aes(y = 0.5*..density..),
                 alpha=0.5,position='identity',binwidth=0.1) +
  ylab("Density") +
  xlab("3' UTR Length") +
  scale_x_continuous(breaks = log(c(10, 100, 1000)),
                     labels = c(10, 100, 1000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = c(0.22, 0.8),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black")) +
  scale_fill_discrete(name = "Fully\nbacked up", labels = c("Yes", "No"))
dev.off()

# Histogram for log-transformed abundance.
fragile.ab.filename <- paste("Scripts/Figures/fragile_abundance_", today.date, ".png", sep = "")
png(fragile.ab.filename, height = 500, width = 500)
ggplot(
  data = ribohits.data,
  aes(
    x = log(Abundance),
    fill = factor(Fragile)
  )
) +
  geom_histogram(aes(y = 0.5*..density..),
                 alpha=0.5,position='identity',binwidth=0.5) +
  ylab("Density") +
  xlab("Abundance") +
  scale_x_continuous(breaks = log(c(0.005, 1, 100, 10000)),
                     labels = c(0, 1, 100, 10000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = c(0.22, 0.8),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black")) +
  scale_fill_discrete(name = "Fully\nbacked up", labels = c("Yes", "No"))
dev.off()

# Histogram for log-transformed in-frame extension lengths.
# First making a data frame with all the extension lengths in one column.
ext.lengths <- data.frame("Fragile" = rep(ribohits.data$Fragile, 3),
                          "Ribohits.Binary" = rep(ribohits.data$Ribohits.Binary, 3),
                          "Length" = c(ribohits.data$Length.0,
                                       ribohits.data$Length.p1,
                                       ribohits.data$Length.p2))
fragile.length.filename <- paste("Scripts/Figures/fragile_length_", today.date, ".png", sep = "")
png(fragile.length.filename, height = 500, width = 500)
ggplot(
  data = ext.lengths,
  aes(
    x = log(Length),
    fill = factor(Fragile)
  )
) +
  geom_histogram(aes(y = 0.5*..density..),
                 alpha=0.5,position='identity',binwidth=0.2) +
  ylab("Density") +
  xlab("Extension length") +
  scale_x_continuous(breaks = log(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  theme_bw(base_size = 28) +
  theme(legend.position = c(0.8025, 0.8375),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black")) +
  scale_fill_discrete(name = "Fully\nbacked up", labels = c("Yes", "No"))
dev.off()
