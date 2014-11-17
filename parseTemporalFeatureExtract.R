# This script will read in the temporal features for 50K gene probes, and parse into one matrix!

setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/timeseries_features")
files = list.files(pattern="*.Rda")

# We will need a matrix of files (genes) by features (160)

rows = c()
load(files[1])
featurematrix = ts$features

for (f in 1:length(files)){
  cat("Processing",f,"of",length(files),"\n")
  load(files[f])
  featurematrix = rbind(featurematrix,ts$features)
  cat("Finishing index",f,"\n")
}

# Save to file
res = list(featurematrix=featurematrix,lastf=37537)
save(res,file="/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/timeseries_matrix_47808_37537.Rda")

