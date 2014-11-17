# This script will visualize results of the developmental enrichment analysis

# Load the result
load("/scratch/users/vsochat/DATA/ALLEN/UberDisorder/developmentRegions1Results/SigDisorderRegionGenes.Rda")

# For each timepoint, save a list of disorders
timepoints = unique(results$TIMEPOINTS)

# The number of disorders at each timepoint
disordercount = c()

# The labels for each timepoint
disorderlabels = c()

for (t in 1:length(timepoints)){
  disorders = unique(results$DISORDER[which(results$TIMEPOINT == timepoints[t])])
  disordercount = c(disordercount,length(disorders))
  disorderlabels = paste(sort(disorders),collapse="\n")
}

# Now plot!
barplot(disordercount)
