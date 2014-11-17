# This script will convert the expression data to timecourses, and calculate temporal features
# We will likely need to run this for each gene - 48K genes by 3702 regions!

library("e1071")

args <- commandArgs(TRUE)
gene = args[1]
outfile = args[2]  # Output file

# Load the timecourses - a list with gene as index (allseries)
load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/genes_47808_timeseries.Rda")

# Define different development periods
tall = colnames(allseries[[1]])
tutero = tall[1:13]
tchild = tall[14:21]
tadul = tall[22:31]
groups = list(all=tall,utero=tutero,child=tchild,adult=tadul)

rowlabels = c()

# We will define the timeseries for each gene, for each region. If we want to average over them at the end, we can do that too
alldata = allseries[[gene]][,groups$all]
interpolated_all = c()

# Now iterate through regions
for (r in 1:nrow(alldata)){
  regiondata = alldata[r,]
  # If we have enough timepoints
  if (length(which(!is.na(regiondata))) >= 8) {
    collabels = names(regiondata)
    # This will hold a row of features for the gene_region across all groups
    label = gsub("[/]|[\\]|[)]|[(]","",gsub(" ","_",paste(gene,rownames(alldata)[r],sep="_")))
    cat("Adding region",label,"\n")
    rowlabels = c(rowlabels,label)
    # Estimate the earliest timepoint, if we don't have it
    if (is.na(regiondata[1])){
      uterotps = regiondata[groups[["utero"]]][!is.na(regiondata[groups[["utero"]]])]
      regiondata[1] = mean(uterotps[length(uterotps)-3:length(uterotps)])   
    }
    # Estimate the latest timepoint, if we don't have it
    if (is.na(regiondata[length(regiondata)])){
      adulttps = regiondata[groups[["adult"]]][!is.na(regiondata[groups[["adult"]]])]
      regiondata[length(regiondata)] = mean(adulttps[length(adulttps)-3:length(adulttps)])   
    }
    # Interpolate missing values
    interpolated = approx(seq(1,length(regiondata)),regiondata,xout=seq(1,length(regiondata)))
   interpolated_all = rbind(interpolated_all,interpolated$y)
  }
}
rownames(interpolated_all) = rowlabels
colnames(interpolated_all) = collabels
save(interpolated_all,file=outfile)
