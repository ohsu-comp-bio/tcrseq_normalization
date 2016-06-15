## Look at lengths of V and J alignments. We want to determine how many of our alignments are good alignments and how
## many may be due to off-target alignment. We're using V alignment length as well as total alignment length to
## determine this. If an alignment length (defined as start of V alignment to end of J alignment) is outside of the
## range of 170-240 basepairs, it is flagged as potential off-target alignment. If the length is outside our range,
## but the V alignment itself is longer than 30 bp, we drop the off-target flag and keep it as a correct alignment.

## objects with clone prefix mean alignments that successfully assembled, defined by having a clone id in the alignment file

##############
### Set up ###
##############

##  Dependencies
library(dplyr)
library(ggplot2)
library(grid)

##  Functions to call later to extract lengths from table
extract.ind <- function(x) {
  start <- as.numeric(strsplit(x, split = '\\|')[[1]][4])
  end <- as.numeric(strsplit(x, split = '\\|')[[1]][5])
  length <- end - start
  row <- c(start, end, length)
  return(row)
} # extract.V


##  Arguments
arguments <- commandArgs(trailingOnly = T)
align.file <- arguments[1]
output.dir <- arguments[2]
##  for testing
align.file <- "~/Desktop/OHSU/tcr_spike/data/equiv_DNA151124LC/export_align/DNA151124LC_1_S1_alignment_exported.txt"

############
### Body ###
############

###
### Data Wrangling and Calculations
###

##  Initialize empty output dataframes for full and subsetted data
summary.df <- NULL
bucket.df <- NULL
clone.summary.df <- NULL
clone.bucket.df <- NULL

##  Read in data
align.data <- read.table(align.file, sep = '\t', header = T,
                        na.strings = c('', ' '), stringsAsFactors = F)

##  Extract sample ID and batch
file.name <- unlist(strsplit(align.file, split = "/"))
file.name <- file.name[length(file.name)]
batch <- unlist(strsplit(file.name, split = "_"))[1]
sample.id <- unlist(strsplit(file.name, split = "_"))[3]

##  Extract V alignment
V.lengths <- t(apply(align.data["Best.V.alignment"], 1, function(x) extract.ind(x)))
colnames(V.lengths) <- c("V.start", "V.end", "V.length")

##  Extract J alignment
J.lengths <- t(apply(align.data["Best.J.alignment"], 1, function(x) extract.ind(x)))
colnames(J.lengths) <- c("J.start", "J.end", "J.length")

##  Combine
align.lengths <- cbind(V.lengths, J.lengths)
align.lengths <- as.data.frame(align.lengths, stringsAsFactors = F)

##  Calculate Total Alignment Length
align.lengths$tot.length <- align.lengths$J.end - align.lengths$V.start

##  Combine to original data
align.data <- cbind(align.data, align.lengths)

##  Subset to include only assembled alignments
clone.data <- align.data[complete.cases(align.data$Clone.Id),]

##  Calculate total alignments (for QC later)
total.alignments <- length(align.data[,1])
total.align.clones <- length(clone.data[,1])

## Begin summary data.frame
summary.df$total.alignments <- total.alignments
clone.summary.df$total.alignments <- total.align.clones

###
###  Buckets of alignment lengths for V and J
###

###   Alignments
##  V
bucket.df$V0.30 <- length(align.data[align.data$V.length <= 30,1])
bucket.df$V31.60 <- length(align.data[align.data$V.length > 30 & align.data$V.length <= 60,1])
bucket.df$V61.90 <- length(align.data[align.data$V.length > 60 & align.data$V.length <= 90,1])
bucket.df$V91.120 <- length(align.data[align.data$V.length > 90 & align.data$V.length <= 120,1])
bucket.df$V121.150 <- length(align.data[align.data$V.length > 120 & align.data$V.length <= 150,1])
bucket.df$V151. <- length(align.data[align.data$V.length > 150,1])

##  J
bucket.df$J0.30 <- length(align.data[align.data$J.length <= 30,1])
bucket.df$J31.60 <- length(align.data[align.data$J.length > 30 & align.data$J.length <= 60,1])
bucket.df$J61.90 <- length(align.data[align.data$J.length > 60 & align.data$J.length <= 90,1])
bucket.df$J91.120 <- length(align.data[align.data$J.length > 90 & align.data$J.length <= 120,1])
bucket.df$J121.150 <- length(align.data[align.data$J.length > 120 & align.data$J.length <= 150,1])
bucket.df$J151. <- length(align.data[align.data$J.length > 150,1])

###   Successfully Assembled Alignments (Clones)
##  V
clone.bucket.df$V0.30 <- length(clone.data[clone.data$V.length <= 30,1])
clone.bucket.df$V31.60 <- length(clone.data[clone.data$V.length > 30 & clone.data$V.length <= 60,1])
clone.bucket.df$V61.90 <- length(clone.data[clone.data$V.length > 60 & clone.data$V.length <= 90,1])
clone.bucket.df$V91.120 <- length(clone.data[clone.data$V.length > 90 & clone.data$V.length <= 120,1])
clone.bucket.df$V121.150 <- length(clone.data[clone.data$V.length > 120 & clone.data$V.length <= 150,1])
clone.bucket.df$V151. <- length(clone.data[clone.data$V.length > 150,1])

##  J
clone.bucket.df$J0.30 <- length(clone.data[clone.data$J.length <= 30,1])
clone.bucket.df$J31.60 <- length(clone.data[clone.data$J.length > 30 & clone.data$J.length <= 60,1])
clone.bucket.df$J61.90 <- length(clone.data[clone.data$J.length > 60 & clone.data$J.length <= 90,1])
clone.bucket.df$J91.120 <- length(clone.data[clone.data$J.length > 90 & clone.data$J.length <= 120,1])
clone.bucket.df$J121.150 <- length(clone.data[clone.data$J.length > 120 & clone.data$J.length <= 150,1])
clone.bucket.df$J151. <- length(clone.data[clone.data$J.length > 150,1])


###
###  Subsetting based on various criteria
###

##  If total alignment length is outside of this range, it is a potentially off-target alignment
align.outside.range <- align.data[align.data$tot.length < 170 | align.data$tot.length > 240, ]
##  Rename
off.target.1 <- align.outside.range
rm(align.outside.range)
##  What percentage of alignments are considered bad at this point?
percent.bad.1 <- round(length(off.target.1[,1]) / total.alignments * 100, digits = 1)

##  If total alignment length is within this range, it is a potentially correct alignment
align.in.range <- align.data[align.data$tot.length >= 170 & align.data$tot.length <= 240,]
##  Rename
correct.1 <- align.in.range
rm(align.in.range)
##  What percentage of alignments are considered good at this point?
percent.good.1 <- round(length(correct.1[,1]) / total.alignments * 100, digits = 1)

##  At this point, we have the most general division between "good" and "bad" alignments. Theoretically though, there are most likely false positives as well as false negatives. Is there something that we can define as a correct alignment to use as a checkpoint here?

## Add to summary
summary.df$percent.outside.of.range <- percent.bad.1
summary.df$percent.in.range <- percent.good.1


##  Now we want to extract from the "bad" alignments any that may actually be "good" alignments and add them to the "good"
##  Any alignment previously flagged as bad, but has an alignment length greater than 30 will be switched to "good".
off.target.2 <- off.target.1[off.target.1$V.length <= 30,]

correct.2 <- rbind(correct.1, off.target.1[off.target.1$V.length > 30,])

##  What percent of total alignments are these?
percent.bad.2 <- round(length(off.target.2[,1]) / total.alignments * 100, digits = 1)
percent.good.2 <- round(length(correct.2[,1]) / total.alignments * 100, digits = 1)

##  Add to summary
summary.df$pct.out.range.short.v <- percent.bad.2
summary.df$pct.in.range.or.long.v <- percent.good.2

##  Combine with buckets

total.df <- data.frame(summary.df, bucket.df)

###################################################################

###  Do the whole thing over using clone data
###
###  Subsetting based on various criteria
###

##  If total alignment length is outside of this range, it is a potentially off-target alignment
clone.outside.range <- clone.data[clone.data$tot.length < 170 | clone.data$tot.length > 240, ]
##  Rename
off.target.1 <- clone.outside.range
rm(clone.outside.range)
##  What percentage of clonements are considered bad at this point?
percent.bad.1 <- round(length(off.target.1[,1]) / total.align.clones * 100, digits = 1)

##  If total alignment length is within this range, it is a potentially correct alignment
clone.in.range <- clone.data[clone.data$tot.length >= 170 & clone.data$tot.length <= 240,]
##  Rename
correct.1 <- clone.in.range
rm(clone.in.range)
##  What percentage of alignments are considered good at this point?
percent.good.1 <- round(length(correct.1[,1]) / total.align.clones * 100, digits = 1)

##  At this point, we have the most general division between "good" and "bad" alignments. Theoretically though, there are most likely false positives as well as false negatives. Is there something that we can define as a correct alignment to use as a checkpoint here?

## Add to summary
clone.summary.df$percent.outside.of.range <- percent.bad.1
clone.summary.df$percent.in.range <- percent.good.1


##  Now we want to extract from the "bad" alignments any that may actually be "good" alignments and add them to the "good"
##  Any alignment previously flagged as bad, but has an alignment length greater than 30 will be switched to "good".
off.target.2 <- off.target.1[off.target.1$V.length <= 30,]

correct.2 <- rbind(correct.1, off.target.1[off.target.1$V.length > 30,])

##  What percent of total alignments are these?
percent.bad.2 <- round(length(off.target.2[,1]) / total.align.clones * 100, digits = 1)
percent.good.2 <- round(length(correct.2[,1]) / total.align.clones * 100, digits = 1)

##  Add to summary
clone.summary.df$pct.out.range.short.v <- percent.bad.2
clone.summary.df$pct.in.range.or.long.v <- percent.good.2

##  Combine with buckets

clone.total.df <- data.frame(clone.summary.df, clone.bucket.df)

###############
### Outputs ###
###############

##  Write summary table
summary.name <- paste(batch, sample.id, "alignment.length.qc", sep = "_")
write.table(total.df,
            file = paste(output.dir, summary.name, sep = ""),
            quote = FALSE,
            sep = ',',
            row.names = FALSE)

##  Write clone summary table
clone.summary.name <- paste(catch, sample.id, "assembled.align.length.qc", sep = "_")
write.table(clone.total.df,
            file = paste(output.dir, clone.summary.name, sep = ''),
            quote = FALSE,
            sep = ',',
            row.names = FALSE)

##  Write good alignments
good.name <- paste(batch, sample.id, "good.alignments.txt", sep = "_")
write.table(correct.2,
            file = paste(output.dir, good.name, sep = ""),
            quote = FALSE,
            sep = '\t',
            row.names = FALSE)
##  Write bad alignments
bad.name <- paste(batch, sample.id, "bad.alignments.txt", sep = "_")
write.table(off.target.2,
            file = paste(output.dir, bad.name, sep = ""),
            quote = FALSE,
            sep = '\t',
            row.names = FALSE)

