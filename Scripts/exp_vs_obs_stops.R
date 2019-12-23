# Calculating expected vs observed probability of stops.

# Import data.
isd.inframe <- read.table("Data/2011ScerSGDCorrectedUTRs.tsv", header = T, sep = "\t", stringsAsFactors = F)
ribohits.data <- read.table("Data/ribohits_data_11-13-19.tab", header = T, sep = "\t", stringsAsFactors = F)

# Calculating nucleotide content to get the probability of having no stops (so we do not bother with C,
# because is never part of a stop codon).
library(tidyverse)

# Simulate via random permutation and return sequences that lack an in-frame backup stop codon.
simulate.nobackup <- function(nt.seq, n = 10000){
  require(stringi)
  require(Biostrings)
  shuffled.seqs <- replicate(n, stri_rand_shuffle(nt.seq))
  dna.set.seqs <- DNAStringSet(shuffled.seqs)
  translated.seqs <- as.character(translate(dna.set.seqs))
  nostop.index <- !grepl("\\*", translated.seqs)
  stopless.seqs <- shuffled.seqs[nostop.index]
  return(stopless.seqs)
}

shift.frame <- function(nt.seq, nt.shift){
  # +1 shift, set nt.shift = 1
  # +2 shift, set nt.shift = 2
  frameshift <- 1 + nt.shift
  shifted.seq <- base::substring(nt.seq, frameshift)
  return(shifted.seq)
}

freq.nostop <- function(nt.seq.set){
  nt.set.dna <- DNAStringSet(nt.seq.set)
  aa.set <- as.character(translate(nt.set.dna))
  nostop.index <- !grepl("\\*", aa.set)
  nt.set.nostop <- nt.seq.set[nostop.index]
  prob.nostop <- length(nt.set.nostop) / length(nt.seq.set)
  return(prob.nostop)
}

prob.nostop.given.nostop <- function(nt.seq, n, nt.shift){
  stopless.seqs <- simulate.nobackup(nt.seq, n)
  stopless.shifted <- shift.frame(stopless.seqs, nt.shift)
  stopless.prob.nostop <- freq.nostop(stopless.shifted)
  return(stopless.prob.nostop)
}

# Sampling from the probabilities to see if they're a success (lack backup) via simulation
# or failure (have backup), like independent bernoulli trials. I'll then count the successes and see
# if they total the observed. I'll then repeat this to see how often they total the observed.
nobackup.df <- tibble(
  "Gene" = ribohits.data$Gene,
  "Length.0" = ribohits.data$Length.0,
  "Length.p2" = ribohits.data$Length.p2,
  "Length.p1" = ribohits.data$Length.p1,
  "Prob.0.p2" = rep(NA, length(ribohits.data$Gene)),
  "Prob.0.p1" = rep(NA, length(ribohits.data$Gene))
)
nobackup.df
set.seed(255)
for (i in 1:length(nobackup.df$Gene)) {
  if (is.na(nobackup.df$Length.0[i])) {
    nobackup.df$Prob.0.p2[i] <- prob.nostop.given.nostop(isd.inframe$UTR[i], 100000, 2)
  }
}
sum(nobackup.df$Prob.0.p2, na.rm = T)
test.sample <- rep(NA, length(nobackup.df[is.na(nobackup.df$Length.0),]$Gene))
for (i in 1:length(nobackup.df[is.na(nobackup.df$Length.0),]$Gene)) {
  test.sample[i] <- rbinom(1, 1, prob = nobackup.df[is.na(nobackup.df$Length.0),]$Prob.0.p2[i])
}
test.sample
sum(test.sample)

bernoulli.sample <- function(probs, n = 1, size = 1){
  # Takes a vector of probabilities and treats each as an independent Bernoulli random variable.
  # Runs a single trial of every Bernoulli random variable, and returns a vector of 1s and 0s,
  # ordered respective to the input probabilities, of success and failure.
  success.vector <- rep(NA, length(probs))
  for (i in 1:length(probs)) {
    success.vector[i] <- rbinom(n = n, size = size, prob = probs[i])
  }
  return(success.vector)
}

sum(bernoulli.sample(nobackup.df[is.na(nobackup.df$Length.0),]$Prob.0.p2))
length(bernoulli.sample(nobackup.df[is.na(nobackup.df$Length.0),]$Prob.0.p2))

length(ribohits.data[is.na(ribohits.data$Length.0) & is.na(ribohits.data$Length.p2), "Gene"])

# Doing that Bernoulli sample 1000 times, seeing the proportion of values that are at least 80.
exp.success.0.p2 <- replicate(10000, sum(bernoulli.sample(nobackup.df[is.na(nobackup.df$Length.0),]$Prob.0.p2)))
exp.success.0.p2[exp.success.0.p2 >= 80]

# Making a function to any frame.
calc.prob.nostop <- function(nt.seq.set, length.set, n, nt.shift){
  # Uses NAs from length to find the sequences without backup stop codons.
  stopifnot(length(nt.seq.set) == length(length.set))
  prob.vector <- rep(NA, length(nt.seq.set))
  for (i in 1:length(nt.seq.set)) {
    if (is.na(length.set[i])) {
      prob.vector[i] <- prob.nostop.given.nostop(nt.seq.set[i], n, nt.shift)
    }
  }
  return(prob.vector)
}
test.prob.vect <- calc.prob.nostop(isd.inframe$UTR, ribohits.data$Length.0, 10, 1)
test.prob.vect
sum(test.prob.vect, na.rm = T)

# Finding the other probabilities. Note that I add "A" or "AA" to switch to the correct starting
# frame. This simulates a +1/+2 shift at a TAA stop codon (e.g., adding "A" is like a +2 frame
# shift, where the last nucleotide of the stop codon becomes the first nucleotide of the first
# codon of the extension). However, it doesn't matter which stop codon we pick for this exercise,
# because the first nucleotide of the extension of a frameshift can never be a stop codon (that
# is, any codon starting in A, G, AA, AG, or GA is never a stop codon). Since this script is ONLY
# calculating whether a stop codon is present or not, the only thing that matters is whether a
# codon is a stop codon or not. So, we can get away with just using the TAA start codon, because
# all that matters is that we gaurantee that the first codon is not a stop.
set.seed(255)
nobackup.df$Prob.0.p2 <- calc.prob.nostop(
  isd.inframe$UTR, ribohits.data$Length.0, 100000, 2
)
set.seed(255)
nobackup.df$Prob.0.p1 <- calc.prob.nostop(
  isd.inframe$UTR, ribohits.data$Length.0, 100000, 1
)
nobackup.df$Prob.p2.0 <- calc.prob.nostop(
  paste0("A", isd.inframe$UTR, sep = ""), ribohits.data$Length.p2, 100000, 1
)
nobackup.df$Prob.p2.p1 <- calc.prob.nostop(
  paste0("A", isd.inframe$UTR, sep = ""), ribohits.data$Length.p2, 100000, 2
)
nobackup.df$Prob.p1.0 <- calc.prob.nostop(
  paste0("AA", isd.inframe$UTR, sep = ""), ribohits.data$Length.p1, 100000, 2
)
nobackup.df$Prob.p1.p2 <- calc.prob.nostop(
  paste0("AA", isd.inframe$UTR, sep = ""), ribohits.data$Length.p1, 100000, 1
)

# Testing proportion of samples with at least the same number of missing second backups as observed.
set.seed(255)
simulations.n <- 10000
obs.nobackup.0.p2 <- length(ribohits.data[is.na(ribohits.data$Length.0) & is.na(ribohits.data$Length.p2), "Gene"])
exp.success.0.p2 <-
  replicate(simulations.n, sum(bernoulli.sample(nobackup.df[is.na(nobackup.df$Length.0),]$Prob.0.p2)))
length(exp.success.0.p2[exp.success.0.p2 >= obs.nobackup.0.p2]) / length(exp.success.0.p2)
exp.success.p2.0 <-
  replicate(simulations.n, sum(bernoulli.sample(nobackup.df[is.na(nobackup.df$Length.p2),]$Prob.p2.0)))
length(exp.success.p2.0[exp.success.p2.0 >= obs.nobackup.0.p2]) / length(exp.success.p2.0)

obs.nobackup.0.p1 <- length(ribohits.data[is.na(ribohits.data$Length.0) & is.na(ribohits.data$Length.p1), "Gene"])
exp.success.0.p1 <-
  replicate(simulations.n, sum(bernoulli.sample(nobackup.df[is.na(nobackup.df$Length.0),]$Prob.0.p1)))
length(exp.success.0.p1[exp.success.0.p1 >= obs.nobackup.0.p1]) / length(exp.success.0.p1)
exp.success.p1.0 <-
  replicate(simulations.n, sum(bernoulli.sample(nobackup.df[is.na(nobackup.df$Length.p1),]$Prob.p1.0)))
length(exp.success.p1.0[exp.success.p1.0 >= obs.nobackup.0.p1]) / length(exp.success.p1.0)

obs.nobackup.p2.p1 <- length(ribohits.data[is.na(ribohits.data$Length.p2) & is.na(ribohits.data$Length.p1), "Gene"])
exp.success.p2.p1 <-
  replicate(simulations.n, sum(bernoulli.sample(nobackup.df[is.na(nobackup.df$Length.p2),]$Prob.p2.p1)))
length(exp.success.p2.p1[exp.success.p2.p1 >= obs.nobackup.p2.p1]) / length(exp.success.p2.p1)
exp.success.p1.p2 <-
  replicate(simulations.n, sum(bernoulli.sample(nobackup.df[is.na(nobackup.df$Length.p1),]$Prob.p1.p2)))
length(exp.success.p1.p2[exp.success.p1.p2 >= obs.nobackup.p2.p1]) / length(exp.success.p1.p2)