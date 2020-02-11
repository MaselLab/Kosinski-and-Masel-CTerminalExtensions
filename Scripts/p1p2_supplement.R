# Supplement frameshifts.

# Load packages.
library(tidyverse)

# Reading data.
ribohits.data <- read_tsv(file = "Data/ribohits_data_11-13-19.tab")
ribohits.data

# Today's date.
today.date <- "1-23-20"
working.dir <- "Scripts/Figures/"
med.ribo <- log10(median(ribohits.data[ribohits.data$Fragile == 0,]$Ribohits.summed + 0.5, na.rm = T))
no.ribo <- log10(0.5)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Excluding genes lacking at least one backup stop codon.
ribo.nofragile <- ribohits.data[ribohits.data$Fragile == 0,]

# Plus 1 extensions first.
p1.orig.lm <- lm(
  data = ribo.nofragile,
  formula = sqrt(ISD.p1.iupred2) ~ log10(Ribohits.summed + 0.5) + log10(Length.p1),
  weights = Length.p1
)
summary(p1.orig.lm)
p1.orig.lrt <- drop1(p1.orig.lm, test = "Chisq")
p1.orig.lrt$`Pr(>Chi)`

p1.rsa.lm <- lm(
  data = ribo.nofragile,
  formula = RSA.p1 ~ log10(Ribohits.summed + 0.5) + log10(Length.p1),
  weights = Length.p1
)
summary(p1.rsa.lm)
p1.rsa.lrt <- drop1(p1.rsa.lm, test = "Chisq")
p1.rsa.lrt$`Pr(>Chi)`

p1.disprop.lm <- lm(
  data = ribo.nofragile,
  formula = Disprop.p1 ~ log10(Ribohits.summed + 0.5) + log10(Length.p1),
  weights = Length.p1
)
summary(p1.disprop.lm)
p1.disprop.lrt <- drop1(p1.disprop.lm, test = "Chisq")
p1.disprop.lrt$`Pr(>Chi)`

# R2 values.
p1.nolength.orig.lm <- lm(
  data = ribo.nofragile,
  formula = sqrt(ISD.p1.iupred2) ~ log10(Ribohits.summed + 0.5),
  weights = Length.p1
)
summary(p1.nolength.orig.lm)

p1.nolength.rsa.lm <- lm(
  data = ribo.nofragile,
  formula = RSA.p1 ~ log10(Ribohits.summed + 0.5),
  weights = Length.p1
)
summary(p1.nolength.rsa.lm)

p1.unweighted.rsa.lm <- lm(
  data = ribo.nofragile,
  formula = RSA.p1 ~ log10(Ribohits.summed + 0.5) + log10(Length.p1)
)
summary(p1.unweighted.rsa.lm)

p1.nolength.disprop.lm <- lm(
  data = ribo.nofragile,
  formula = Disprop.p1 ~ log10(Ribohits.summed + 0.5),
  weights = Length.p1
)
summary(p1.nolength.disprop.lm)

p1.unweighted.disprop.lm <- lm(
  data = ribo.nofragile,
  formula = Disprop.p1 ~ log10(Ribohits.summed + 0.5) + log10(Length.p1)
)
summary(p1.unweighted.disprop.lm)

# Plotting.
p1.ribo.isd.plot.filename <- paste(working.dir, "ISD_ribo_p1_", today.date, ".png", sep = "")
png(filename = p1.ribo.isd.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = sqrt(ISD.p1.iupred2),
    weight = Length.p1,
    size = Length.p1
  )
) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("ISD") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = sqrt(c(0, 0.1, 0.25, 0.5, 1)),
                     labels = c("0.00", "0.10", "0.25", "0.50", "1.00")) +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

p1.length.isd.plot.filename <- paste(working.dir, "ISD_length_p1_", today.date, ".png", sep = "")
png(filename = p1.length.isd.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Length.p1),
    y = sqrt(ISD.p1.iupred2),
    color = Ribohits
  )
) +
  stat_function(fun = function(x)p1.orig.lm$coefficients[1]+p1.orig.lm$coefficients[3]*x+no.ribo*p1.orig.lm$coefficients[2],
                geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)p1.orig.lm$coefficients[1]+p1.orig.lm$coefficients[3]*x+med.ribo*p1.orig.lm$coefficients[2],
                geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("ISD") +
  scale_x_continuous(breaks = log10(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  scale_y_continuous(breaks = sqrt(c(0.04, 0.16, 0.36)),
                     labels = c(0.04, 0.16, 0.36),
                     limits = c(0.125, 0.62)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), labels = c("Median", "None"), values = cbPalette[c(2,4)]) +
  theme(legend.position = c(0.75, 0.8),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black"))
dev.off()

p1.ribo.rsa.plot.filename <- paste(working.dir, "RSA_ribo_p1_", today.date, ".png", sep = "")
png(filename = p1.ribo.rsa.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = RSA.p1,
    weight = Length.p1,
    size = Length.p1
  )
) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("RSA") +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

p1.length.rsa.plot.filename <- paste(working.dir, "RSA_length_p1_", today.date, ".png", sep = "")
png(filename = p1.length.rsa.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Length.p1),
    y = RSA.p1,
    color = Ribohits
  )
) +
  stat_function(fun = function(x)p1.rsa.lm$coefficients[1]+p1.rsa.lm$coefficients[3]*x+no.ribo*p1.rsa.lm$coefficients[2],
                geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)p1.rsa.lm$coefficients[1]+p1.rsa.lm$coefficients[3]*x+med.ribo*p1.rsa.lm$coefficients[2],
                geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("RSA") +
  scale_x_continuous(breaks = log10(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), labels = c("Median", "None"), values = cbPalette[c(2,4)]) +
  theme(legend.position = c(0.75, 0.8),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black"))
dev.off()

p1.ribo.disprop.plot.filename <- paste(working.dir, "Disprop_ribo_p1_", today.date, ".png", sep = "")
png(filename = p1.ribo.disprop.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = Disprop.p1,
    weight = Length.p1,
    size = Length.p1
  )
) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("Disorder propensity") +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

p1.length.disprop.plot.filename <- paste(working.dir, "Disprop_length_p1_", today.date, ".png", sep = "")
png(filename = p1.length.disprop.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Length.p1),
    y = Disprop.p1,
    color = Ribohits
  )
) +
  stat_function(fun = function(x)p1.disprop.lm$coefficients[1]+p1.disprop.lm$coefficients[3]*x+no.ribo*p1.disprop.lm$coefficients[2],
                geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)p1.disprop.lm$coefficients[1]+p1.disprop.lm$coefficients[3]*x+med.ribo*p1.disprop.lm$coefficients[2],
                geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("Disorder propensity") +
  scale_x_continuous(breaks = log10(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), labels = c("Median", "None"), values = cbPalette[c(2,4)]) +
  theme(legend.position = "none")
dev.off()

### Unweighted versions of the plots.
p1.unweighted.isd.plot.filename <- paste(working.dir, "ISD_unweighted_ribo_p1_", today.date, ".png", sep = "")
png(filename = p1.unweighted.isd.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = sqrt(ISD.p1.iupred2)
  )
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
dev.off()

p1.unweighted.rsa.plot.filename <- paste(working.dir, "RSA_unweighted_ribo_p1_", today.date, ".png", sep = "")
png(filename = p1.unweighted.rsa.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = RSA.p1)
) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("RSA") +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

p1.unweighted.disprop.plot.filename <- paste(working.dir, "Disprop_unweighted_ribo_p1_", today.date, ".png", sep = "")
png(filename = p1.unweighted.disprop.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = Disprop.p1)
) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("Disorder propensity") +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

#############################################################################
#############################################################################
# And now plus 2 extensions.
p2.orig.lm <- lm(
  data = ribo.nofragile,
  formula = sqrt(ISD.p2.iupred2) ~ log10(Ribohits.summed + 0.5) + log10(Length.p2),
  weights = Length.p2
)
summary(p2.orig.lm)
p2.orig.lrt <- drop1(p2.orig.lm, test = "Chisq")
p2.orig.lrt$`Pr(>Chi)`

p2.grafted.lm <- lm(
  data = ribo.nofragile,
  formula = sqrt(ISD.p2.grafted) ~ log10(Ribohits.summed + 0.5) + log10(Length.p2),
  weights = Length.p2
)
summary(p2.grafted.lm)
p2.grafted.lrt <- drop1(p2.grafted.lm, test = "Chisq")
p2.grafted.lrt$`Pr(>Chi)`

p2.rsa.lm <- lm(
  data = ribo.nofragile,
  formula = RSA.p2 ~ log10(Ribohits.summed + 0.5) + log10(Length.p2),
  weights = Length.p2
)
summary(p2.rsa.lm)
p2.rsa.lrt <- drop1(p2.rsa.lm, test = "Chisq")
p2.rsa.lrt$`Pr(>Chi)`

p2.disprop.lm <- lm(
  data = ribo.nofragile,
  formula = Disprop.p2 ~ log10(Ribohits.summed + 0.5) + log10(Length.p2),
  weights = Length.p2
)
summary(p2.disprop.lm)
p2.disprop.lrt <- drop1(p2.disprop.lm, test = "Chisq")
p2.disprop.lrt$`Pr(>Chi)`

# R2 values.
p2.nolength.orig.lm <- lm(
  data = ribo.nofragile,
  formula = sqrt(ISD.p2.iupred2) ~ log10(Ribohits.summed + 0.5),
  weights = Length.p2
)
summary(p2.nolength.orig.lm)

p2.nolength.grafted.lm <- lm(
  data = ribo.nofragile,
  formula = sqrt(ISD.p2.grafted) ~ log10(Ribohits.summed + 0.5),
  weights = Length.p2
)
summary(p2.nolength.grafted.lm)

p2.nolength.rsa.lm <- lm(
  data = ribo.nofragile,
  formula = RSA.p2 ~ log10(Ribohits.summed + 0.5),
  weights = Length.p2
)
summary(p2.nolength.rsa.lm)

p2.nolength.disprop.lm <- lm(
  data = ribo.nofragile,
  formula = Disprop.p2 ~ log10(Ribohits.summed + 0.5),
  weights = Length.p2
)
summary(p2.nolength.disprop.lm)
################################
p2.nolength.unweighted.orig.lm <- lm(
  data = ribo.nofragile,
  formula = sqrt(ISD.p2.iupred2) ~ log10(Ribohits.summed + 0.5)
)
summary(p2.nolength.unweighted.orig.lm)

p2.nolength.unweighted.grafted.lm <- lm(
  data = ribo.nofragile,
  formula = sqrt(ISD.p2.grafted) ~ log10(Ribohits.summed + 0.5)
)
summary(p2.nolength.unweighted.grafted.lm)

p2.unweighted.rsa.lm <- lm(
  data = ribo.nofragile,
  formula = RSA.p2 ~ log10(Ribohits.summed + 0.5) + log10(Length.p2)
)
summary(p2.unweighted.rsa.lm)

p2.nolength.unweighted.rsa.lm <- lm(
  data = ribo.nofragile,
  formula = RSA.p2 ~ log10(Ribohits.summed + 0.5)
)
summary(p2.nolength.unweighted.rsa.lm)

p2.unweighted.disprop.lm <- lm(
  data = ribo.nofragile,
  formula = Disprop.p2 ~ log10(Ribohits.summed + 0.5) + log10(Length.p2)
)
summary(p2.unweighted.disprop.lm)

p2.nolength.unweighted.disprop.lm <- lm(
  data = ribo.nofragile,
  formula = Disprop.p2 ~ log10(Ribohits.summed + 0.5)
)
summary(p2.nolength.unweighted.disprop.lm)

# Plotting.
p2.ribo.isd.plot.filename <- paste(working.dir, "ISD_ribo_p2_", today.date, ".png", sep = "")
png(filename = p2.ribo.isd.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = sqrt(ISD.p2.iupred2),
    weight = Length.p2,
    size = Length.p2
  )
) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("ISD") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = sqrt(c(0, 0.1, 0.25, 0.5, 1)),
                     labels = c("0.00", "0.10", "0.25", "0.50", "1.00")) +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

p2.length.isd.plot.filename <- paste(working.dir, "ISD_length_p2_", today.date, ".png", sep = "")
png(filename = p2.length.isd.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Length.p2),
    y = sqrt(ISD.p2.iupred2),
    color = Ribohits
  )
) +
  stat_function(fun = function(x)p2.orig.lm$coefficients[1]+p2.orig.lm$coefficients[3]*x+no.ribo*p2.orig.lm$coefficients[2],
                geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)p2.orig.lm$coefficients[1]+p2.orig.lm$coefficients[3]*x+med.ribo*p2.orig.lm$coefficients[2],
                geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("ISD") +
  scale_x_continuous(breaks = log10(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  scale_y_continuous(breaks = sqrt(c(0.04, 0.16, 0.36)),
                     labels = c(0.04, 0.16, 0.36),
                     limits = c(0.125, 0.65)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), labels = c("Median", "None"), values = cbPalette[c(2,4)]) +
  theme(legend.position = c(0.75, 0.8),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black"))
dev.off()

p2.ribo.isd.graftedplot.filename <- paste(working.dir, "ISD_ribo_p2_g_", today.date, ".png", sep = "")
png(filename = p2.ribo.isd.graftedplot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = sqrt(ISD.p2.grafted),
    weight = Length.p2,
    size = Length.p2
  )
) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("ISD") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = sqrt(c(0, 0.1, 0.25, 0.5, 1)),
                     labels = c("0.00", "0.10", "0.25", "0.50", "1.00")) +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

p2.length.isd.graftedplot.filename <- paste(working.dir, "ISD_length_p2_g_", today.date, ".png", sep = "")
png(filename = p2.length.isd.graftedplot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Length.p2),
    y = sqrt(ISD.p2.grafted),
    color = Ribohits
  )
) +
  stat_function(fun = function(x)p2.grafted.lm$coefficients[1]+p2.grafted.lm$coefficients[3]*x+no.ribo*p2.grafted.lm$coefficients[2],
                geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)p2.grafted.lm$coefficients[1]+p2.grafted.lm$coefficients[3]*x+med.ribo*p2.grafted.lm$coefficients[2],
                geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("ISD") +
  scale_x_continuous(breaks = log10(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  scale_y_continuous(breaks = sqrt(c(0.04, 0.16, 0.36)),
                     labels = c(0.04, 0.16, 0.36),
                     limits = c(0.125, 0.65)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), labels = c("Median", "None"), values = cbPalette[c(2,4)]) +
  theme(legend.position = "none")
dev.off()

p2.ribo.rsa.plot.filename <- paste(working.dir, "RSA_ribo_p2_", today.date, ".png", sep = "")
png(filename = p2.ribo.rsa.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = RSA.p2,
    weight = Length.p2,
    size = Length.p2
  )
) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("RSA") +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

p2.length.rsa.plot.filename <- paste(working.dir, "RSA_length_p2_", today.date, ".png", sep = "")
png(filename = p2.length.rsa.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Length.p2),
    y = RSA.p2,
    color = Ribohits
  )
) +
  stat_function(fun = function(x)p2.rsa.lm$coefficients[1]+p2.rsa.lm$coefficients[3]*x+no.ribo*p2.rsa.lm$coefficients[2],
                geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)p2.rsa.lm$coefficients[1]+p2.rsa.lm$coefficients[3]*x+med.ribo*p2.rsa.lm$coefficients[2],
                geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("RSA") +
  scale_x_continuous(breaks = log10(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), labels = c("Median", "None"), values = cbPalette[c(2,4)]) +
  theme(legend.position = c(0.75, 0.8),
        legend.background = element_rect(fill = "white", size = 1, linetype = "solid", color = "black"))
dev.off()

p2.ribo.disprop.plot.filename <- paste(working.dir, "Disprop_ribo_p2_", today.date, ".png", sep = "")
png(filename = p2.ribo.disprop.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = Disprop.p2,
    weight = Length.p2,
    size = Length.p2
  )
) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("Disorder propensity") +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none")
dev.off()

p2.length.disprop.plot.filename <- paste(working.dir, "Disprop_length_p2_", today.date, ".png", sep = "")
png(filename = p2.length.disprop.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Length.p2),
    y = Disprop.p2,
    color = Ribohits
  )
) +
  stat_function(fun = function(x)p2.disprop.lm$coefficients[1]+p2.disprop.lm$coefficients[3]*x+no.ribo*p2.disprop.lm$coefficients[2],
                geom = "line", aes(color = "No"), size = 2) +
  stat_function(fun = function(x)p2.disprop.lm$coefficients[1]+p2.disprop.lm$coefficients[3]*x+med.ribo*p2.disprop.lm$coefficients[2],
                geom = "line", aes(color = "Yes"), size = 2) +
  xlab("Extension length (AAs)") +
  ylab("Disorder propensity") +
  scale_x_continuous(breaks = log10(c(1, 10, 100)),
                     labels = c(1, 10, 100)) +
  theme_bw(base_size = 28) +
  scale_color_manual(limits = c("Yes", "No"), labels = c("Median", "None"), values = cbPalette[c(2,4)]) +
  theme(legend.position = "none")
dev.off()

### Unweighted version of the +2 plots.
p2.ribo.isd.plot.filename <- paste(working.dir, "ISD_unweighted_ribo_p2_", today.date, ".png", sep = "")
png(filename = p2.ribo.isd.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = sqrt(ISD.p2.iupred2))
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
dev.off()

p2.ribo.isd.graftedplot.filename <- paste(working.dir, "ISD_unweighted_ribo_p2_g_", today.date, ".png", sep = "")
png(filename = p2.ribo.isd.graftedplot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = sqrt(ISD.p2.grafted))
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
dev.off()

p2.ribo.rsa.plot.filename <- paste(working.dir, "RSA_unweighted_ribo_p2_", today.date, ".png", sep = "")
png(filename = p2.ribo.rsa.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = RSA.p2
    )
) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("RSA") +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28)
dev.off()

p2.ribo.disprop.plot.filename <- paste(working.dir, "Disprop_unweighted_ribo_p2_", today.date, ".png", sep = "")
png(filename = p2.ribo.disprop.plot.filename, width = 6, height = 6, units = "in", res = 350)
ggplot(
  data = ribo.nofragile,
  aes(
    x = log10(Ribohits.summed + 0.5),
    y = Disprop.p2
  )
) +
  geom_point(size = 0.5) +
  geom_smooth(method = "loess") +
  geom_smooth(method = "lm", color = "red") +
  xlab("Ribohits") +
  ylab("Disorder propensity") +
  scale_x_continuous(breaks = log10(c(0.5, 1.5, 10.5, 100.5, 1000.5)),
                     labels = c(0, 1, 10, 100, 1000)) +
  theme_bw(base_size = 28)
dev.off()