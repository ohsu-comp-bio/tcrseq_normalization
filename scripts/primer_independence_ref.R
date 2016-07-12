
# This is the external R script that is referenced by primer_independence_regression.Rmd

###
### Read in Spike Count data, combine, melt, normalize
###

### Input: Directory of 25-bp spike count files
### Output: list of  1) original data frame prior to modification
###                  2) melted dataframe containing 1 row for each VJ combination 
###                     for each sample and the corresponding spike count. Also contains
###                     one column for each VJ combo and a 1 if that row has same combo,
###                     a 0, if not.
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
  
  # Collapse data frame to 1 count column
  melt.count.df <- melt(count.df, id.vars = c("V", "J"))
  
  # Add pseudo-variable of V/J combos
  melt.count.df$combos <- paste(melt.count.df$V, melt.count.df$J, sep = '')
  
  # Add a column for each variable combo
  combos <- paste(count.df$V, count.df$J, sep = '')
  for (i in 1:length(combos)){
    curr.combo <- vector(mode = "numeric", length = length(melt.count.df[,1]))
    curr.name <- combos[i]
    melt.count.df[,curr.name] <- curr.combo
  } # for i
  
  # Mark as 1 every row where the combos value equals the row.name
  for (i in 6:length(melt.count.df[1,])){
    col.combo <- colnames(melt.count.df)[i]
    for (j in 1:length(melt.count.df[,1])){
      if (melt.count.df[j,5] == col.combo){
        melt.count.df[j,i] <- 1
      } # if
    } # for j
  } # for i
  
  # Take log2 of count values due to geometric distribution
  log2.melt.count.df <- melt.count.df
  log2.melt.count.df$value <- log2(melt.count.df$value + 1)
  
  # Change colnames for step compatibility
  colnames(log2.melt.count.df) <- gsub("-", "_", colnames(log2.melt.count.df))
  colnames(melt.count.df) <- gsub("-", "_", colnames(melt.count.df))
  
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
  log2.wo.int.coeff <- summary(log2.wo.int.lm)$coefficients
  print(paste("Log2 Without Interaction R^2:", log2.wo.int.adj.r2, sep = ' '))
  # With interaction
  log2.with.int.lm <- lm(value ~ V * J, batch$log2)
  log2.with.int.adj.r2 <- round(summary(log2.with.int.lm)$adj.r.squared, digits = 4)
  log2.with.int.coeff <- summary(log2.with.int.lm)$coefficients
  print(paste("Log2 With Interaction R^2:", log2.with.int.adj.r2, sep = ' '))
  
  
  return(list("log2.wo.int" = log2.wo.int.lm, "log2.wo.int.r2" = log2.wo.int.adj.r2, "log2.wo.int.coeff" = log2.wo.int.coeff,
              "log2.with.int" = log2.with.int.lm, "log2.with.int.r2" = log2.with.int.adj.r2, "log2.with.int.coeff" =log2.with.int.coeff))
} # make.models(batch)

###
### Step function to determine best model
###

### Input: 1) List of data frames that is output by read.data()
###        2) Vector containing all the column names of VJ combinations to be tested
### Output: List containing the model and it's summary that was chosen by the step function

steps <- function(batch, vjcombos, k.val){

    # Null
    null <- lm(value ~ 1, data = batch$log2)
    
    
    # Full
    full.mod <- paste("value ~ V + J +", paste(vjcombos, collapse = " + "), sep = '')
    full.mod <- as.formula(full.mod)
    full <- lm(full.mod, data = batch$log2)
    
    # Step forward and backward
    print("Both: ")
    both <- step(null, scope = list(upper = full), data = batch$log2,
                 direction = "both", k = k.val)
    both.summ <- summary(both)
    
    return(list("both" = both, "both.summary" = both.summ))
} # steps(batch)

# steps.w.subset <- function(batch, subset){
#   null <- lm(value ~ 2, data = batch$log2)
#   full <- lm(value ~ V + J + V5-J2-5, data = batch$log2)
#              
#   both <- step(full, scope = list(lower = null, upper = full),
#              direction = "both")
#   both.summ <- summary(both)
#   return(list("both" = both, "both.summary" = both.summ))
# } # steps.w.subset
	


###
### Find VJ combinations that significantly influence model
###

### Input: 1) list of dataframes from read.data()
###        3) significance value to use

### Output:

# vj.combos.to.use <- function(lm.coeff, sig.value = numeric()){
#   
#   if (!missing(sig.value)){
#     # Extract from coefficient table based on significance avlue
#     vj.combos <- which(lm.coeff[,4] <= sig.value)
#     vj.combos <- lm.coeff[vj.combos,]
#     # Get just names and remove Intercept
#     vj.combos <- rownames(vj.combos)
#     vj.combos <- vj.combos[2:length(vj.combos)]
#     # Change to match column names from data frame
#     vj.combos <- gsub("^V|:J", '', vj.combos)
#     vj.combos <- gsub("-", "_", vj.combos)
#     vj.combos <- gsub("_$", "-", vj.combos)
#     return(vj.combos)
#   } else{
#     vj.combos <- rownames(lm.coeff)
#     vj.combos <- vj.combos[2:length(vj.combos)]
#     # Change to match column names from data frame
#     vj.combos <- gsub("^V|:J", '', vj.combos)
#     vj.combos <- gsub("-", "_", vj.combos)
#     vj.combos <- gsub("_$", "-", vj.combos)
#     
#   } # if
#   
# }# vj.combos.to.use
  
  # Extract coefficients with p-vaues above specified threshold value
  #model.coeff <- step.results[[2]]$coefficients
  #insig.p.values <- which(model.coeff[,4] > sig.value)
  #sig.p.values <- which(model.coeff[,4] <= sig.value)
  
  #new.coeff <- model.coeff[insig.p.values,]
  #removed.coeff <- model.coeff[sig.p.values,]
  
  # Get variable IDs
  #new.coeff.ids <- gsub("combos", '', row.names(new.coeff))
  
  # Subset count data
  #new.melt.data <- batch$log2[which(batch$log2$combos %in% new.coeff.ids),]
  
  # Return
  #return(list("log2" = new.melt.data, "removed combos" = removed.coeff))
#} # remove.sig.variables(batch, step.results, sig.value)
