# Need to combine all of the qc data into a single data frame.


##  Arguments

arguments <- commandArgs(trailingOnly = T)
qc.dir <- arguments[1]
clone.qc.dir <- arguments[2]

##  List files
qc.files <- list.files(qc.dir)
clone.qc.files <- list.files(clone.qc.dir)

##  Order files
qc.files <- qc.files[order(as.numeric(gsub(".*_S|_align.*", '', qc.files)))]
clone.qc.files <- clone.qc.files[order(as.numeric(gsub(".*_S|_assembl.*", '', clone.qc.files)))]

## Metadata
batch <- unlist(strsplit(qc.files[1], split = "_"))[1]

##  Empty data frame
total.qc <- data.frame(sample.id=character(), total.alignments=numeric(), percent.outside.of.range=numeric(),
	    				      percent.in.range=numeric(), pct.out.range.short.v=numeric(),
					      pct.in.range.or.long.v=numeric(), V0.30 = numeric(), V31.60 = numeric(), V61.90 = numeric(),
					      V91.120 = numeric(), V121.150 = numeric(), V151. = numeric(), J0.30 = numeric(), J31.60 = numeric(), 
					      J61.90 = numeric(), J91.120 = numeric(), J121.150 = numeric(), J151. = numeric())

clone.total.qc <- total.qc

##  Combine files

for (i in 1:length(qc.files)){
    curr.qc <- read.csv(paste(qc.dir, qc.files[i], sep = ''), stringsAsFactors = F, header = T)

    sample.id <- unlist(strsplit(qc.files[i], split = "_"))[2]

    row <- data.frame(sample.id, curr.qc[1,])
#    print(row)
    total.qc <- rbind(total.qc, row, stringsAsFactors = F)

}  # for 

for (i in 1:length(clone.qc.files)){
  curr.clone.qc <- read.csv(paste(clone.qc.dir, clone.qc.files[i], sep = ''), stringsAsFactors = F, header = T)
  
  sample.id <- unlist(strsplit(clone.qc.files[i], split = "_"))[2]
  
  row <- data.frame(sample.id, curr.clone.qc[1,])
  #    print(row)
  clone.total.qc <- rbind(clone.total.qc, row, stringsAsFactors = F)
  
}  # for 


##  Write output

output.name <- paste(batch, "aggregate.align.length.qc.txt", sep = "_")
clone.output.name <- paste(batch, "aggregate.clone.align.length.qc.txt", sep = "_")

write.table(total.qc,
		file = paste(qc.dir, output.name, sep = ''),
		sep = '\t',
		quote = F,
		row.names = F)

write.table(clone.total.qc,
            file = paste(clone.qc.dir, clone.output.name, sep = ''),
            sep = '\t',
            quote = F,
            row.names = F)