# New Figure 2 (histogram)

# Load ribohits data.
ribohits.data <- read.table("Data/ribohits_data_11-13-19.tab", header = T, sep = "\t", stringsAsFactors = F)

# Load packages.
library(tidyverse)

# Global variables.
today.date <- "1-24-20"

# Histogram of 3' UTR ribohits coverage.
histogram.ribo.filename <- paste("Scripts/Figures/histogram_ribohits_", today.date, ".png", sep = "")
png(histogram.ribo.filename, width = 15, height = 7, units = "in", res = 350)
ggplot(
  data = ribohits.data,
  aes(
    x = log10(Ribohits.summed + 0.5),
    fill = as.character(Fragile)
  )
) +
  geom_histogram(breaks = c(seq(-0.301030, 3.7, by = 0.5)), alpha = 0.7) +
  #geom_smooth() +
  xlab("3' UTR Ribohits") +
  scale_x_continuous(breaks = c(seq((-0.301030 + 0.19897)/2, 3.7, by = 0.5)),
                     labels = c("0", "1-4", "5-15", "16-49", "50-157", "158-499", "500-1580", "1581-2199")) +
  #scale_y_continuous(breaks = sqrt(c(0, 0.1, 0.25, 0.5, 1)),
  #                   labels = c("0.00", "0.10", "0.25", "0.50", "1.00")) +
  theme_bw(base_size = 28) +
  theme(legend.position = c(0.75, 0.8),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black")) +
  scale_fill_discrete(name = "Fully\nbacked up", labels = c("Yes", "No"))
dev.off()
