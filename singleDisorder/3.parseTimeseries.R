# This script will read in the interpolated timeseries, and parse
setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/interpolated_timeseries")
files = list.files(pattern="*_timeseries.Rda")

# We will parse rows into one matrix
timeseries = c()
genes = c()

for (f in 1:length(files)){
  cat("Processing",f,"of",length(files),"\n")
  gene_name = strsplit(files[f],"_timeseries.Rda")[[1]][1]
  load(files[f])
  # Here is the number of genes
  num = dim(interpolated_all)[1]
  genes = c(genes,rep(gene_name,num))
  timeseries = rbind(timeseries,interpolated_all)
}

# Save to file
interpolated = list(timeseries=timeseries,genes=genes)
save(interpolated,file="/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/interpolated_timeseries.Rda")



