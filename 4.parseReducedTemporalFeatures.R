# This script will read in the reduced temporal features for ~50K gene probes, and parse into one matrix!  We hope that this matrix is much smaller than the original with ~700K rows

setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/geneRegionSimilarity")
files = list.files(pattern="*_reduced.Rda")

# We will parse rows into one matrix
temporalmatrix = c()
rows = c()

# We should also save a complete list to lookup regions in each
regions = list()

for (f in 1:length(files)){
  cat("Processing",f,"of",length(files),"\n")
  load(files[f])
  temporalmatrix = rbind(temporalmatrix,result$reduced)
  rows = c(rows,names(result$regions))
  regions = c(regions,result$regions)
}

# Add rownames
rownames(temporalmatrix) = rows

# Save to file
temporal = list(matrix=temporalmatrix,rowlookup=regions,rows=rows)
save(temporal,file="/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/reduced_feature_matrix_51641.Rda")

