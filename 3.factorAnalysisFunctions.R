##############################################
# Functions for the analysis of questionnaires
##############################################


# function to select the biggest value of the factor loading matrix by row
likefadiagram <- function(facload){
  bin.load <- matrix(0,dim(facload)[1],dim(facload)[2])
  for(i in 1:dim(facload)[1]){
    bin.load[i,which.max(abs(facload[i,]))] <- 1
  }
  bin.load.NA <- replace(bin.load,bin.load==0,NA)
  return(list(binary=bin.load, load=facload * bin.load, load.NA=facload * bin.load.NA))
}

# function to select values larger than percentile 'p'
likefadiagram.p <- function(facload,p){
  bin.load <- matrix(0,dim(facload)[1],dim(facload)[2])
  for(i in 1:dim(facload)[1]){
    cte <- quantile(abs(facload[i,]),p)
    bin.load[i,which(abs(facload[i,])>cte)] <- 1
  }
  bin.load.NA <- replace(bin.load,bin.load==0,NA)
  return(list(binary=bin.load,load=facload*bin.load, load.NA=facload * bin.load.NA))
}

# function to select values larger than percentile 'p' (global)
likefadiagram.per <- function(facload,p){
  cte <- quantile(abs(facload),p)
  bin.load <- matrix(0,dim(facload)[1],dim(facload)[2])
  for(i in 1:dim(facload)[1]){
    bin.load[i,which(abs(facload[i,])>cte)] <- 1
  }
  bin.load.NA <- replace(bin.load,bin.load==0,NA)
  return(list(binary=bin.load,load=facload*bin.load, load.NA=facload * bin.load.NA))
}
