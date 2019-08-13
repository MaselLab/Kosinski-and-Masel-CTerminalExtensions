# Table 1 script and other non-model tests.

# Load ribohits data.
ribohits.data <- read.table("Data/ribohits_data_tai_8-12-19.tsv", header = T, sep = "\t", stringsAsFactors = F)

# Generate table.
fragile.ribo <- table(ribohits.data[, c("Ribohits.Binary", "Fragile")])
fragile.ribo

# Chi-square test on table.
fragile.chi <- chisq.test(fragile.ribo)
fragile.chi$p.value

# Fragile vs. non-fragile protein t-tests for abundance, 3' UTR length,
# extension length, and C-terminal ISD.
ab.fragile <- 
  t.test(log(ribohits.data[ribohits.data$Fragile == 1, "Abundance"]),
         log(ribohits.data[ribohits.data$Fragile == 0, "Abundance"]),
         var.equal = T)
ab.fragile$p.value
fragile.ab.mean <- exp(mean(log(ribohits.data[ribohits.data$Fragile == 1, "Abundance"]), na.rm = T))
nonfragile.ab.mean <- exp(mean(log(ribohits.data[ribohits.data$Fragile == 0, "Abundance"]), na.rm = T))
fragile.ab.fold.change <- nonfragile.ab.mean / fragile.ab.mean
fragile.ab.fold.change

utr.fragile <-
  t.test(log(ribohits.data[ribohits.data$Fragile == 1,]$Length.3UTR),
         log(ribohits.data[ribohits.data$Fragile == 0,]$Length.3UTR),
         var.equal = T)
utr.fragile$p.value

in.frame.fragile <-
  t.test(log(ribohits.data[ribohits.data$Fragile == 1,]$Length.0),
         log(ribohits.data[ribohits.data$Fragile == 0,]$Length.0),
         var.equal = T)
in.frame.fragile

m1.fragile <-
  t.test(log(ribohits.data[ribohits.data$Fragile == 1,]$Length.m1),
         log(ribohits.data[ribohits.data$Fragile == 0,]$Length.m1),
         var.equal = T)
m1.fragile

m2.fragile <-
  t.test(log(ribohits.data[ribohits.data$Fragile == 1,]$Length.m2),
         log(ribohits.data[ribohits.data$Fragile == 0,]$Length.m2),
         var.equal = T)
m2.fragile

cterm.isd.fragile <-
  t.test(sqrt(ribohits.data[ribohits.data$Fragile == 1,]$ISD.Last10.iupred2),
         sqrt(ribohits.data[ribohits.data$Fragile == 0,]$ISD.Last10.iupred2),
         var.equal = T)
cterm.isd.fragile

# Showing that in-frame extensions are shorter.
in.frame.vs.m1.ext <-
  t.test(log(ribohits.data$Length.0),
         log(ribohits.data$Length.m1),
         var.equal = T)
in.frame.vs.m1.ext$p.value

in.frame.vs.m2.ext <-
  t.test(log(ribohits.data$Length.0),
         log(ribohits.data$Length.m2),
         var.equal = T)
in.frame.vs.m2.ext$p.value

# Showing that in-frame extensions are shorter for non-fragile and fragile proteins.
# Non-fragile.
ext.0.m1.nofragile <- 
  with(ribohits.data[ribohits.data$Fragile == 0,],
       t.test(log(Length.0), log(Length.m1), var.equal = T))
ext.0.m2.nofragile <- 
  with(ribohits.data[ribohits.data$Fragile == 0,],
       t.test(log(Length.0), log(Length.m2), var.equal = T))
ext.0.m1.nofragile$p.value
ext.0.m2.nofragile$p.value
# Fragile.
ext.0.m1.fragile <- 
  with(ribohits.data[ribohits.data$Fragile == 1,],
       t.test(log(Length.0), log(Length.m1), var.equal = T))
ext.0.m2.fragile <- 
  with(ribohits.data[ribohits.data$Fragile == 1,],
       t.test(log(Length.0), log(Length.m2), var.equal = T))
ext.0.m1.fragile$p.value
ext.0.m2.fragile$p.value

# Are C-termini of genes with detectable ribohits more disordered than genes without detectable ribohits?
ribo.vs.isd.last10 <-
  t.test(sqrt(ribohits.data[ribohits.data$Ribohits.Binary == 1 & ribohits.data$Fragile == 0,]$ISD.Last10.iupred2),
         sqrt(ribohits.data[ribohits.data$Ribohits.Binary == 0 & ribohits.data$Fragile == 0,]$ISD.Last10.iupred2),
         var.equal = T)
ribo.vs.isd.last10

# Is geometric mean tAI (for the first 4 codon positions in in-frame extensions) higher in non-fragile than
# in fragile genes?
t.test(ribohits.data[ribohits.data$Fragile == 0,]$tai.4.log,
       ribohits.data[ribohits.data$Fragile == 1,]$tai.4.log,
       var.equal = T)
