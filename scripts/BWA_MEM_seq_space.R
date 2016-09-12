# We are using the BWA-MEM aligner to check our MiXCR results. One step in the process requires that we make a subsetted sequence
# space that only contains the TCR-Beta genes (instead of the entire mouse genome). We can do this by using custom fasta files
# that contain only the sequence for V, D, or J TCR-Beta genes as input and then iteratively pasting together different combinations.
# These VDJ combinations that we create will not have the junctions that give normal T-Cells their diversity (and allow us to count
# clonotypes), but we believe the aligner will reduce all of the alignment scores equally because none of the reference sequences
# will have them and all of our sequences will.

### Inputs
# Have three input files: 
#       V, D, and J fasta files. 
#       Header line contains >chr:start-end

### Outputs
# 1. Fasta file of #V * #D * #J sequences and a header line identifying which V, D, and J created it
# 2. Tab-separated file containing the V-end, J-start, and total insert length, for each VDJ

##############
### Set up ###
##############

### Read in fasta files
v.reads <- readFasta("~/BWA-MEM/seq_space_construction/fasta/trbv.fa")
d.reads <- readFasta("~/BWA-MEM/seq_space_construction/fasta/trbd.fa")
d.reads <- d.reads[1]
j.reads <- readFasta("~/BWA-MEM/seq_space_construction/fasta/trbj.fa")

#################
### Calculate ###
#################

## Create empty variables
output <- NULL
insert.lengths <- NULL

## Iterate through each of the V genes, and extract sequence, ID, start, and end
for (i in 1:length(v.reads)){
  curr.v.seq <- v.reads@sread[i]
  curr.v.id <- v.reads@id[i]
  curr.v.end <- as.numeric(gsub(".*-|_.*", '', curr.v.id))
  curr.v.start <- as.numeric(gsub(".*:|-.*", '', curr.v.id))
  ## During each V gene iteration, iterate through each D gene
  ## extract sequence and ID
  ## paste D to end of V
  for (j in 1:length(d.reads)){
    curr.d.seq <- d.reads@sread[j]
    curr.d.id <- d.reads@id[j]
    curr.first.join <- paste(curr.v.seq, curr.d.seq, sep = '')
    ## During each VD iteration, iterate through each J gene
    ## extract sequence, ID, start, and end
    ## paste J to end of D
    for(k in 1:length(j.reads)){
      curr.j.seq <- j.reads@sread[k]
      curr.j.id <- j.reads@id[k]
      curr.j.start <-  as.numeric(gsub(".*:|-.*", '', curr.j.id))
      curr.j.end <- as.numeric(gsub('[0-9]{1}:[0-9]{5,10}-|_.*', '', curr.j.id))
      curr.final.join <- paste(curr.first.join, curr.j.seq, sep = '')

      ## Turn into a DNA sequence object and add custom header
      curr.final.seq <- DNAStringSet(x = curr.final.join)
      curr.final.seq@ranges@NAMES <- paste("6:", curr.v.start, "-", curr.j.end, "_V:", i, "_D:", j, "_J:", k, sep = '')

      ## Write to output file
      writeFasta(curr.final.seq, "~/BWA-MEM/ref/VDJ.one.copy.fa", mode = "a")

      ## Calculate insert length - end of V gene to beginning of J gene
      ## Append to data frame
      curr.ins.length <- curr.j.start - curr.v.end
      ins.row <- c(curr.v.end, curr.j.start, curr.ins.length)
      insert.lengths <- rbind(insert.lengths, ins.row)

      ## Below is for iteration checking
      ## cat(i, j, k, '\n')
    } # for k
  } # for j
} # for i

## Change to data frame and add column names
insert.lengths <- data.frame(insert.lengths)
colnames(insert.lengths) <- c("V.end", "J.start", "Insert.Length")

## Write output
write.table(insert.lengths, file = "~/BWA-MEM/ref/insert.lengths.txt", quote = F, row.names = F, col.names = T,
            sep = '\t')

## Read fasta back in again to check against reads
# check <- readFasta("~/BWA-MEM/ref/new.seq.space.fa")
# output[224,1] == check@sread[112]
# output[2,1] == check@sread[1]
# output[100,1] == check@sread[50]
# output[576,1] == check@sread[288]

## Calculate insert length summary statistics
ins.mean <- mean(insert.lengths$Insert.Length)
ins.stdev <- sd(insert.lengths$Insert.Length)
ins.min <- min(insert.lengths$Insert.Length)
ins.max <- max(insert.lengths$Insert.Length)

ins.length.summary <- c("mean" = ins.mean, "stdev" = ins.stdev, "min" = ins.min, "max" = ins.max)

