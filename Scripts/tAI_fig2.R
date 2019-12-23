# tAI analysis

# Load packages.
library(tidyverse)

# Load ribohits and sequence data.
seq.data <- read.table("Data/2011_Scer_SGD_UTRs_7-8-19.txt", header = T, sep = "\t", stringsAsFactors = F)
ribohits.data <- read.table("Data/ribohits_data_11-13-19.tab", header = T, sep = "\t", stringsAsFactors = F)
head(ribohits.data)
head(seq.data)

# Load tAI data.
tai.data <- read.table("Scripts/tAI/gene_tAI.txt", header = T, stringsAsFactors = F, sep = ",")
head(tai.data)

# Importing the necessary information from the ribohits data (length and fragile).
identical(ribohits.data$Gene, tai.data$Header)
tai.data$Fragile <- ribohits.data$Fragile
tai.data$Length.0 <- ribohits.data$Length.0
identical(ribohits.data$Gene, seq.data$Gene)
seq.data$ISD.0 <- ribohits.data$ISD.0.iupred2

# Function for calculating geometric mean tAI up to each position.
tai.position.function <- function(
  start.pos = 1,
  tai.df = tai.data,
  seq.df = seq.data,
  frame.name = "In-frame",
  ext.aa.col = "Ext_AA",
  length.col = "Length.0",
  check.if.na.col = "ISD.0"
) {
  if (is.na(seq.df[1, check.if.na.col]) == F) {
    if (str_length(seq.df[1, ext.aa.col]) + 1 - start.pos > 1){
      position.df <- data.frame("tAI" = 
                                  stack(tai.df[1,
                                               (str_length(seq.df$ORF_AA[1]) + start.pos + 2):
                                                 (str_length(seq.df$ORF_AA[1]) +
                                                    str_length(seq.df[1, ext.aa.col]) + 2)]
                                  )[,1],
                                "Position" = start.pos:str_length(seq.df[1, ext.aa.col]),
                                "Gene" = rep(tai.df$Header[1], str_length(seq.df[1, ext.aa.col]) + 1 - start.pos),
                                "Fragile" = rep(tai.df$Fragile[1], str_length(seq.df[1, ext.aa.col]) + 1 - start.pos),
                                #"Ribohits" = rep(tai.df$Ribohits[1], str_length(seq.df[1, ext.aa.col]) + 1 - start.pos),
                                "Length" = rep(tai.df[1, length.col], str_length(seq.df[1, ext.aa.col]) + 1 - start.pos),
                                "Frame" = rep(frame.name, str_length(seq.df[1, ext.aa.col]) + 1 - start.pos),
                                row.names = NULL
      )
    } else if (str_length(seq.df[1, ext.aa.col]) + 1 - start.pos == 1) {
      position.df <- data.frame(
        "tAI" = tai.df[1, (str_length(seq.df$ORF_AA[1]) + start.pos + 2):
                         (str_length(seq.df$ORF_AA[1]) + str_length(seq.df[1, ext.aa.col]) + 2)],
        "Position" = start.pos,
        "Gene" = tai.df$Header[1],
        "Fragile" = tai.df$Fragile[1],
        #"Ribohits" = tai.df$Ribohits[1],
        "Length" = tai.df[1, length.col],
        "Frame" = frame.name,
        row.names = NULL
      )
    } else if (str_length(seq.df[1, ext.aa.col]) + 1 - start.pos < 1) {
      position.df <- data.frame(
        "tAI" = NA,
        "Position" = start.pos,
        "Gene" = tai.df$Header[1],
        "Fragile" = tai.df$Fragile[1],
        #"Ribohits" = tai.df$Ribohits[1],
        "Length" = tai.df[1, length.col],
        "Frame" = frame.name,
        row.names = NULL
      )
    } 
  } else if (is.na(seq.df[1, check.if.na.col]) == T) {
    position.df <- data.frame(
      "tAI" = NA,
      "Position" = start.pos,
      "Gene" = tai.df$Header[1],
      "Fragile" = tai.df$Fragile[1],
      #"Ribohits" = tai.df$Ribohits[1],
      "Length" = tai.df[1, length.col],
      "Frame" = frame.name,
      row.names = NULL
    )
  } else {
    print("What happened?")
    print(tai.df$Header[i])
  }
  for (i in 2:length(tai.df$Header)) {
    if (is.na(seq.df[i, check.if.na.col]) == F) {
      if (str_length(seq.df[i, ext.aa.col]) + 1 - start.pos > 1){
        temp.df <- data.frame(
          "tAI" = stack(
            tai.df[i,
                   (str_length(seq.df$ORF_AA[i]) + start.pos + 2):
                     (str_length(seq.df$ORF_AA[i]) + str_length(seq.df[i, ext.aa.col]) + 2)]
          )[,1],
          "Position" = start.pos:str_length(seq.df[i, ext.aa.col]),
          "Gene" = rep(tai.df$Header[i], str_length(seq.df[i, ext.aa.col]) + 1 - start.pos),
          "Fragile" = rep(tai.df$Fragile[i], str_length(seq.df[i, ext.aa.col]) + 1 - start.pos),
          #"Ribohits" = rep(tai.df$Ribohits[i], str_length(seq.df[i, ext.aa.col]) + 1 - start.pos),
          "Length" = rep(tai.df[i, length.col], str_length(seq.df[i, ext.aa.col]) + 1 - start.pos),
          "Frame" = rep(frame.name, str_length(seq.df[i, ext.aa.col]) + 1 - start.pos),
          row.names = NULL
        )
        #print(position.df)
        #print(temp.df)
        position.df <- rbind(position.df, temp.df)
      } else if (str_length(seq.df[i, ext.aa.col]) + 1 - start.pos == 1) {
        temp.df <- data.frame(
          "tAI" = tai.df[i, (str_length(seq.df$ORF_AA[i]) + start.pos + 2):
                           (str_length(seq.df$ORF_AA[i]) + str_length(seq.df[i, ext.aa.col]) + 2)],
          "Position" = start.pos,
          "Gene" = tai.df$Header[i],
          "Fragile" = tai.df$Fragile[i],
          #"Ribohits" = tai.df$Ribohits[i],
          "Length" = tai.df[i, length.col],
          "Frame" = frame.name,
          row.names = NULL
        )
        position.df <- rbind(position.df, temp.df)
      } else {
        print("What happened?")
        print(tai.df$Header[i])
      }
    }
  }
  return(position.df)
}

tai.position.df <- tai.position.function()

# Checking the data frame to make sure all genes except the ones without extensions are there.
length(ribohits.data[is.na(ribohits.data$ISD.0.iupred2), "Gene"])
4082 - length(unique(tai.position.df$Gene))
identical(ribohits.data[!is.na(ribohits.data$ISD.0.iupred2), "Gene"], as.character(unique(tai.position.df$Gene)))

# Making a cumulative mean tAI across codon positions.
tai.position.df$tAI.cumulative.log.mean <- NA
current.gene <- "blank"
current.frame <- "blank"
for (i in 1:length(tai.position.df$Gene)) {
  if (tai.position.df$Frame[i] != current.frame) {
    current.frame <- tai.position.df$Frame[i]
  }
  if (tai.position.df$Gene[i] != current.gene){
    current.gene <- tai.position.df$Gene[i]
    tAI.vector <- vector()
    tAI.mean.vector <- vector()
    tAI.vector <- tai.position.df[tai.position.df$Gene == current.gene &
                                    tai.position.df$Frame == current.frame,]$tAI
    tAI.mean.vector <- cumsum(log(tAI.vector)) / seq_along(tAI.vector)
    for (j in 1:length(tAI.mean.vector)){
      tai.position.df$tAI.cumulative.log.mean[i + j - 1] <- tAI.mean.vector[j]
    }
  }
}

# Checking significance between fragile and non-fragile proteins.
t.test(
  x = tai.position.df[tai.position.df$Position == 4 &
                        tai.position.df$Frame == "In-frame" &
                        tai.position.df$Fragile == 1,]$tAI.cumulative.log.mean,
  y = tai.position.df[tai.position.df$Position == 4 &
                        tai.position.df$Frame == "In-frame" &
                        tai.position.df$Fragile == 0,]$tAI.cumulative.log.mean,
  var.equal = T
)
