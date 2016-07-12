## Companion script for DNA160609_expected_alignment_length


###
### Extract actual alignment lengths from MiXCR export_align table
###

extract.ind <- function(x) {
  start <- as.numeric(strsplit(x, split = '\\|')[[1]][4])
  end <- as.numeric(strsplit(x, split = '\\|')[[1]][5])
  length <- end - start
  row <- c(start, end, length)
  return(row)
} # extract.V


###
### Calculate V, J, and total alignment lengths and place in data frame
###

get.lengths <- function(align.df) {
  
  ## V
  V.lengths <- t(apply(align.df["Best.V.alignment"], 1, function(x) extract.ind(x)))
  colnames(V.lengths) <- c("V.start", "V.end", "V.length")
  
  ## J
  J.lengths <- t(apply(align.df["Best.J.alignment"], 1, function(x) extract.ind(x)))
  colnames(J.lengths) <- c("J.start", "J.end", "J.length")
  
  ## Combine
  align.lengths <- cbind.data.frame(V.lengths, J.lengths)
  
  ## Calcualte total length
  align.lengths$tot.length <- align.lengths$J.end - align.lengths$V.start
  
  return(align.lengths)
}


###
### Determine Alignment lengths within appropriate range
###

## V
compare.lengths.V <- function(align.df, expected.df) {
  
  # Create empty data frame
  in.range.df <- matrix(nrow = length(align.df[,1]), ncol = length(align.df[1,]))
  colnames(in.range.df) <- colnames(align.df)
  out.range.df <- matrix(nrow = length(align.df[,1]), ncol = length(align.df[1,]))
  colnames(out.range.df) <- colnames(align.df)
  
   for (i in 1:length(align.df[,1])) {
    # Get relevant info
    curr.V <- align.df$Best.V.hit[i]
    curr.length <- align.df$V.length[i]
    # Skip if not one of our V's
    if(!(curr.V %in% vs)){
      next
    }
    # Add to new data frame if V is in our range.
    if((curr.length > expected.df[expected.df == curr.V, 3]) 
       & (curr.length < expected.df[expected.df == curr.V, 4])){
      row <- unlist(align.df[i, ])
      in.range.df[i,] <- row
    } else { # Add to separate data frame if outside our range
      row <- unlist(align.df[i,])
      out.range.df[i,] <- row
    }
  } # for
  in.range.df <- as.data.frame(in.range.df)
  out.range.df <- as.data.frame(out.range.df)
  in.range.df <- in.range.df[complete.cases(in.range.df$Read.id),]
  out.range.df <- out.range.df[complete.cases(out.range.df$Read.id),]
  return(list("good" = in.range.df, "bad" = out.range.df))
} # compare.lengths.V

## J
compare.lengths.J <- function(align.df, expected.df) {
  
  # Create empty data frames
  in.range.df <- matrix(nrow = length(align.df[,1]), ncol = length(align.df[1,]))
  colnames(in.range.df) <- colnames(align.df)
  out.range.df <- matrix(nrow = length(align.df[,1]), ncol = length(align.df[1,]))
  colnames(out.range.df) <- colnames(align.df)
  
  for (i in 1:length(align.df[,1])) {
    # Get relevant info
    curr.J <- align.df$Best.J.hit[i]
    curr.length <- align.df$J.length[i]
    # Skip if not one of our J's
    if(!(curr.J %in% js)){
      next
    }
    # Add to new df if J is in our range
    if((curr.length > expected.df[expected.df == curr.J, 3]) 
       & (curr.length < expected.df[expected.df == curr.J, 4])){
      row <- unlist(align.df[i, ])
      in.range.df[i,] <- row
    } else { # Add to separate df if outside our range
      row <- unlist(align.df[i,])
      out.range.df[i,] <- row
    }
  } # for
  in.range.df <- as.data.frame(in.range.df)
  in.range.df <- in.range.df[complete.cases(in.range.df$Read.id),]
  out.range.df <- as.data.frame(out.range.df)
  out.range.df <- out.range.df[complete.cases(out.range.df$Read.id),]
  return(list("good" = in.range.df, "bad" = out.range.df))
} # compare.lengths.J

## Both
compare.lengths.both <- function(align.df, expected.df){
  
  # Create empty data frame
  in.range.df <- matrix(nrow = length(align.df[,1]), ncol = length(align.df[1,]))
  colnames(in.range.df) <- colnames(align.df)
  out.range.df <- matrix(nrow = length(align.df[,1]), ncol = length(align.df[1,]))
  colnames(out.range.df) <- colnames(align.df)
  
  for (i in 1:length(align.df[,1])) {
    # Get relevant info
    curr.V <- align.df$Best.V.hit[i]
    curr.J <- align.df$Best.J.hit[i]
    curr.J.length <- align.df$J.length[i]
    curr.V.length <- align.df$V.length[i]
    # Skip if not one of our V's or our J's
    if(!(curr.J %in% js) | !(curr.V %in% vs)){
      next
    }
    # Add to new df if V and J are in our range
    if (curr.V.length > expected.df[expected.df == curr.V, 3] & curr.V.length < expected.df[expected.df == curr.V, 4] &
        curr.J.length > expected.df[expected.df == curr.J, 3] & curr.J.length < expected.df[expected.df == curr.J, 4]){
      row <- unlist(align.df[i, ])
      in.range.df[i,] <- row
    } else { # add to separate df if not in our range
      row <- unlist(align.df[i,])
      out.range.df[i,] <- row
    }
  } # for
  in.range.df <- as.data.frame(in.range.df)
  in.range.df <- in.range.df[complete.cases(in.range.df$Read.id),]
  out.range.df <- as.data.frame(out.range.df)
  out.range.df <- out.range.df[complete.cases(out.range.df$Read.id),]
  return(list("good" = in.range.df, "bad" = out.range.df))
} # compare.lengths.both


clean.envir.first.stage <- function() {
  rm(align.outside.range, num.outside.range, align.in.range, num.in.range,
     outside.good.v, num.outside.good.v, outside.bad.v, num.outside.bad.v,
     inside.good.v, num.inside.good.v, inside.bad.v, num.inside.bad.v,
     outside.good.v.assembled, num.outside.good.v.assembled,
     outside.bad.v.assembled, num.outside.bad.v.assembled,
     inside.good.v.assembled, num.inside.good.v.assembled,
     inside.bad.v.assembled, num.inside.bad.v.assembled, envir = globalenv())
} # clean.envir.first.stage