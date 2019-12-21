# Figure 3.

# Load ribohits data.
ribohits.data <- read.table("Data/ribohits_data_11-13-19.tab", header = T, sep = "\t", stringsAsFactors = F)

# Load packages.
library(tidyverse)
library(Hmisc)
library(stringr)

# Import sequence data to calculate lengths.
seq.data <- read.table("Data/2011_Scer_SGD_UTRs_7-8-19.txt", header = T, sep = "\t", stringsAsFactors = F)
seq.data$Length.ORF <- NA
for (i in 1:length(seq.data$Gene)) {
  seq.data$Length.ORF[i] <- stringr::str_length(seq.data$ORF_AA[i])
}
identical(seq.data$Gene, ribohits.data$Gene)
ribohits.data$Length.ORF <- seq.data$Length.ORF
ribohits.data$Length.0 <- ribohits.data$Length.0 - 1

# Gather data into ISD groups.
ribo.isd <- ribohits.data[ribohits.data$Fragile == 0,] %>%
  gather("ISD.ORF.iupred2", "ISD.3UTR.iupred2", "ISD.0.iupred2",
                         key = "SeqType", value = "ISD")
ribo.length <- ribohits.data[ribohits.data$Fragile == 0,] %>%
  gather("Length.ORF", "Length.3UTR", "Length.0",
         key = "SeqType", value = "Length")
identical(ribo.isd$Gene, ribo.length$Gene)
for (i in 1:length(ribo.isd$Gene)) {
  if (ribo.isd$SeqType[i] == "ISD.ORF.iupred2") {
    ribo.isd$SeqType[i] = "ORF"
  } else if (ribo.isd$SeqType[i] == "ISD.3UTR.iupred2") {
    ribo.isd$SeqType[i] = "3' UTR"
  } else if (ribo.isd$SeqType[i] == "ISD.0.iupred2") {
    ribo.isd$SeqType[i] = "Ext"
  } else {
    print(ribo.isd$Gene[i])
  }
}
for (i in 1:length(ribo.length$Gene)) {
  if (ribo.length$SeqType[i] == "Length.ORF") {
    ribo.length$SeqType[i] = "ORF"
  } else if (ribo.length$SeqType[i] == "Length.3UTR") {
    ribo.length$SeqType[i] = "3' UTR"
  } else if (ribo.length$SeqType[i] == "Length.0") {
    ribo.length$SeqType[i] = "Ext"
  } else {
    print(ribo.length$Gene[i])
  }
}
identical(ribo.isd$SeqType, ribo.length$SeqType)
ribo.isd$Length <- ribo.length$Length

isd.seqtype <- group_by(ribo.isd, SeqType)
isd.data <- isd.seqtype %>%
  dplyr::summarise(ISD.mean = wtd.mean(sqrt(ISD), weights = Length),
            ISD.sd = sqrt(wtd.var(ISD, weights = Length)),
            ISD.N = length(ISD[!is.na(ISD)]),
            ISD.se = ISD.sd / sqrt(ISD.N))
isd.data$SeqType <- factor(isd.data$SeqType, levels = c("ORF", "3' UTR", "Ext"))
isd.seqtype$SeqType <- factor(isd.seqtype$SeqType, levels = c("ORF", "3' UTR", "Ext"))

# Making figure 2 part A with means and standard errors.
png(filename = "Scripts/Figures/ISD_ORF_UTR_Ext_12-18-19.png", width = 6, height = 6, units = "in", res = 350)
ggplot(data = isd.data,
       aes(x = SeqType,
           y = ISD.mean ^ 2)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = (ISD.mean - ISD.se)^2, ymax = (ISD.mean + ISD.se)^2), width = 0.4, size = 2) +
  scale_x_discrete(breaks = c("ORF", "3' UTR", "Ext"),
                   labels = c("ORF", "3' UTR", "Ext")) +
  xlab("") +
  ylab("ISD") +
  theme_bw(base_size = 28)
dev.off()

# Making figure 2 part A with boxplot.
boxplot.quantiles <- function(x){
  qntls <- quantile(x, probs = c(0.09, 0.25, 0.5, 0.75, 0.91))
  names(qntls) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(qntls)
}

png(filename = "Scripts/Figures/ISD_ORF_UTR_Ext_boxplot_12-18-19.png", height = 400, width = 400)
ggplot(data = isd.seqtype,
       aes(x = SeqType,
           y = sqrt(ISD))) +
  stat_summary(fun.data = boxplot.quantiles, geom = "boxplot") +
  scale_x_discrete(breaks = c("ORF", "3' UTR", "Ext"),
                   labels = c("ORF", "3' UTR", "Ext")) +
  scale_y_continuous(breaks = sqrt(c(0.05, 0.15, 0.30, 0.5)),
                     labels = c("0.05", "0.15", "0.30", "0.50")) +
  xlab("") +
  ylab("ISD") +
  theme_bw(base_size = 28)
dev.off()

# Making figure 2 part b.
abundance.lm <- lm(data = ribohits.data[ribohits.data$Fragile == 0,],
                   formula = sqrt(ISD.Last10.iupred2) ~ log(Abundance))
ab.lm.test <- drop1(abundance.lm, test = "Chisq")
ab.lm.test$`Pr(>Chi)`
ribo.nofragile <- ribohits.data[ribohits.data$Fragile == 0,]
abundance.last10.pred <-
  with(ribo.nofragile,
     predict(abundance.lm, newdata = data.frame("Abundance" = Abundance,
                                                "ISD.Last10.iupred2" = ISD.Last10.iupred2,
                                                "Length.3UTR" = median(Length.3UTR, na.rm = T),
                                                "ISD.grafted.iupred2" = median(ISD.grafted.iupred2, na.rm = T),
                                                "Length.0" = median(Length.0, na.rm = T))))
ribo.nofragile$ab.fit <- abundance.last10.pred

png(filename = "Scripts/Figures/Abundance_ISDLast10_12-18-19.png", width = 6, height = 6, units = "in", res = 350)
ggplot(data = ribo.nofragile,
       aes(y = sqrt(ISD.Last10.iupred2),
           x = log(Abundance))) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", color = "red", se = F, size = 2) +
  ylab("C-terminal ISD") +
  xlab("Protein abundance (ppm)") +
  scale_x_continuous(breaks = log(c(0.007, 1, 100, 10000)),
                     labels = c(0, 1, 100, 10000)) +
  scale_y_continuous(breaks = sqrt(c(0, 0.25, 0.5, 0.75, 1)),
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

# Checking effect size.
predict(abundance.lm, newdata = data.frame("Abundance" = quantile(ribo.nofragile$Abundance, na.rm = T))) ^ 2
