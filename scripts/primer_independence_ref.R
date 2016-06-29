
# This is the external R script that is referenced by primer_independence_regression.Rmd

###
### Populate data frame with appropriate information
###

populate.vj.df <- function(df.of.counts, vj.df){
     for (i in 1:length(df.of.counts[,3])){
          v <- as.character(df.of.counts[i,1])
       	  j <- as.character(df.of.counts[i,2])
          sum <- df.of.counts[i,3]
	  rownum <- which(v == row.names(vj.df))
	  colnum <- which(j == colnames(vj.df))
	  vj.df[rownum,colnum] <- sum
     }  # for
     return(vj.df)
}  # populate.vj.df(df.of.counts, vj.df)

###
### Read in Spike Count data, combine, melt, normalize
###

### Input: Directory of 25-bp spike count files
### Output: list of  1) original data frame prior to modification
###                  2) melted dataframe containing 1 row for each VJ combination 
###                     for each sample and the corresponding spike count.
###                  3) The above data frame, with spike counts normalized using log2

read.data <- function(count.dir){
  
  # Read in files and sort by sample number
  all.counts <- list.files(count.dir)
  all.counts <- all.counts[order(as.numeric(gsub(".*_S|.assembled.*", '',
                                                 all.counts)))]
  
  # Read in first file to start aggregate data frame
  count.df <- read.table(file.path(count.dir, all.counts[1]), sep =',', header = T)
  count.df <- count.df[,3:5]
  
  # Comine spike counts for all files into 1 data frame
  # Columns are samples
  # Rows are spikes
  for (i in 2:length(all.counts)){
    curr.df <- read.table(file.path(count.dir, all.counts[i]), sep = ',', 
                          header = T)
    count.df <- cbind(count.df, curr.df$spike.count)
  }   #   for i in 2:length(new.counts)
  colnames(count.df) <- c("V", "J", seq(1:length(all.counts)))
  
  #head(count.df[,1:10], n = 20)
  
  # Collapse data frame to 1 count column
  melt.count.df <- melt(count.df, id.vars = c("V", "J"))
  
  # Add pseudo-variable of V/J combos
  melt.count.df$combos <- paste(melt.count.df$V, melt.count.df$J, sep = '')
  
  # Take log2 of count values due to geometric distribution
  log2.melt.count.df <- melt.count.df
  log2.melt.count.df$value <- log2(melt.count.df$value + 1)
  
  # Divide all values by upper quartile as another normalization method
  #norm.melt.count.df <- melt.count.df
  #norm.melt.count.df$value <- norm.melt.count.df$value /
  #  summary(norm.melt.count.df$value)[5]
  
  # V is all V segments, repeated for each J for each sample
  # J is all J segments, repeated same as V
  # Variable corresponds to sample number
  # Value is log2 of count
  
  return(list("original" = count.df, "melt" = melt.count.df, 
              "log2" = log2.melt.count.df)) #, "third.q" = norm.melt.count.df))
} # read.data(count.df)


###
### Construct models
###

### Input: List of data frames that is output by read.data()
### Output: List of models and R2 values. First is without interaction, and second is with

make.models <- function(batch){
  ### Using Log2 Data
  # Without interaction
  log2.wo.int.lm <- lm(value ~ V + J, batch$log2)
  log2.wo.int.adj.r2 <- round(summary(log2.wo.int.lm)$adj.r.squared, digits = 4)
  print(paste("Log2 Without Interaction R^2:", log2.wo.int.adj.r2, sep = ' '))
  # With interaction
  log2.with.int.lm <- lm(value ~ V * J, batch$log2)
  log2.with.int.adj.r2 <- round(summary(log2.with.int.lm)$adj.r.squared, digits = 4)
  print(paste("Log2 With Interaction R^2:", log2.with.int.adj.r2, sep = ' '))
  
  
  ### Using Third Quartile data
  # Additive
  #third.q.add.lm <- lm(value ~ V + J, batch$third.q)
  #third.q.add.adj.r2 <- round(summary(third.q.add.lm)$adj.r.squared, digits = 4)
  #print(paste("Third Quartile Additive R^2:", third.q.add.adj.r2, sep = ' '))
  
  # Multiplicative
  #third.q.mult.lm <- lm(value ~ V * J, batch$third.q)
  #third.q.mult.adj.r2 <- round(summary(third.q.mult.lm)$adj.r.squared, digits = 4)
  #print(paste("Third Quartile Multiplicative R^2:", third.q.mult.adj.r2, sep = ' '))
  
  
  return(list("log2.wo.int" = log2.wo.int.lm, "log2.wo.int.r2" = log2.wo.int.adj.r2,
              "log2.with.int" = log2.with.int.lm, "log2.with.int.r2" = log2.with.int.adj.r2)) #,
              #"third.q.add" = third.q.add.lm, "third.q.add.r2" = third.q.add.adj.r2,
              #"third.q.mult" = third.q.add.lm, "third.q.mult.r2" = third.q.mult.adj.r2))
} # make.models(batch)

###
### Step function to determine best model
###

### Input: List of data frames that is output by read.data() or find.sig.variables()
### Output: List containing the model and it's summary that was chosen by the step function

steps <- function(batch){
#  per.norm.method <- function(batch.method){
#    print(deparse(substitute(batch.method)))
    # Null
    null <- lm(value ~ 1, data = batch$log2)
    # Full
    full <- lm(value ~ V + J + combos, data = batch$log2)
#    full <- lm(value  ~ V + J + V:J, data = batch$log2)
    # Step forward and backward
    print("Both: ")
    both <- step(full, scope = list(lower = null, upper = full),
                 direction = "both")
    both.summ <- summary(both)
    
 #   return(list("both" = both, "both.summary" = both.summ))
 # } # per.norm.method(batch.method)
  
#  batch.log2 <- per.norm.method(batch$log2)
  #batch.third.q <- per.norm.method(batch$third.q)
  
#  return(batch.log2)
    return(list("both" = both, "both.summary" = both.summ))
  #return(list("log2" = batch.log2,
  #            "third.quartile" = batch.third.q))
} # steps(batch)
	


###
### Extract non-significant variables from model and remove from dataset
###

### Input: 1) List of data frames output by read.table()
###        2) List of model and summary output by steps()
### Output: List containing newly subsetted dataframe containing only 
###         significant VJ combos and a vector of those VJ combos' identities

remove.unsig.variables <- function(batch, step.results){
  # Extract coefficients with significant P values from model
  model.coeff <- step.results[[2]]$coefficients
  sig.p.values <- which(model.coeff[,4] < 0.05)
  insig.p.values <- which(model.coeff[,4] >= 0.05)
  new.coeff <- model.coeff[sig.p.values,]
  bad.coeff <- model.coeff[insig.p.values,]
  
  # Get variable IDs
  new.coeff.ids <- gsub('combos', '', row.names(new.coeff))
  #bad.coeff.ids <- gsub('combos', '', row.names(bad.coeff))
  
  # Use IDs to subset count data
  new.melt.data <- batch$log2[which(batch$log2$combos %in% new.coeff.ids),]
  
  # Return 
  return(list("log2" = new.melt.data, "sig.ids" = new.coeff.ids, "unsig.coeff" = bad.coeff))
  
} # remove.unsig.variables(batch, step.results)


remove.sig.variables <- function(batch, step.results, sig.value = numeric()){
  
  # Extract coefficients with p-vaues above specified threshold value
  model.coeff <- step.results[[2]]$coefficients
  insig.p.values <- which(model.coeff[,4] > sig.value)
  sig.p.values <- which(model.coeff[,4] <= sig.value)
  
  new.coeff <- model.coeff[insig.p.values,]
  removed.coeff <- model.coeff[sig.p.values,]
  
  # Get variable IDs
  new.coeff.ids <- gsub("combos", '', row.names(new.coeff))
  
  # Subset count data
  new.melt.data <- batch$log2[which(batch$log2$combos %in% new.coeff.ids),]
  
  # Return
  return(list("log2" = new.melt.data, "removed combos" = removed.coeff))
} # remove.sig.variables(batch, step.results, sig.value)

iterate.sig.combos <- function(batch, step.results) {
  
  # Extract coefficients with p-values euqal to < 2e-16
  model.coeff <- step.results[[2]]$coefficients
  p.vals.to.iterate <- model.coeff[which(model.coeff[,4] < 2 * 10 ** -16),]
  
  # Get IDs
  ids.to.iterate <- gsub("combos", '', row.names(p.vals.to.iterate))
  ids.to.iterate <- ids.to.iterate[2:length(ids.to.iterate)]
  
  # Make empty data frame
  combo.remove.compare <- matrix(nrow = length(ids.to.iterate), ncol = 2)
  colnames(combo.remove.compare) <- c("Primer.Combo", "R2")
  
  # Iterate over list, removing each VJ combination, and re-running the step function
  for (i in 1:length(ids.to.iterate)) {
    new.melt.data <- batch$log2[which(!(batch$log2$combos %in% ids.to.iterate[i])),]
    curr.lm <- with(new.melt.data, lm(value ~ combos))
    curr.lm.summary <- summary(curr.lm)
    row <- c(ids.to.iterate[i], as.numeric(curr.lm.summary$adj.r.squared))
    combo.remove.compare[i,] <- row
  } # for
  combo.remove.compare <- as.data.frame(combo.remove.compare, stringsAsFactors = F)
  return(combo.remove.compare)
} # iterate.sig.combos