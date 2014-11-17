# This script will convert the expression data to timecourses, and calculate temporal features
# We will likely need to run this for each gene - 48K genes by 3702 regions!

library("e1071")

args <- commandArgs(TRUE)
gene = args[1]
imgdir = args[2]  # Directory to save GENE_REGION.png images
outfile = args[3]  # Output file

# Load the timecourses - a list with gene as index (allseries)
load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/genes_47808_timeseries.Rda")

# Define different development periods
tall = colnames(allseries[[1]])
tutero = tall[1:13]
tchild = tall[14:21]
tadul = tall[22:31]
groups = list(all=tall,utero=tutero,child=tchild,adult=tadul)

# USER FUNCTIONS ----------------------------------------------------------------------------
shannon.entropy <- function(p) {
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

# TEMPORAL FEATURE EXTRACTION --------------------------------------------------------------
temporal = c() # The matrix of features, gene-region combos in rows, features in columns
rowlabels = c()  # the gene-region combinations
collabels = c() # feature labels

# We will define the timeseries for each gene, for each region. If we want to average over them at the end, we can do that too

alldata = allseries[[gene]][,groups$all]
# Now iterate through regions
for (r in 1:nrow(alldata)){
  regiondata = alldata[r,]
  # If we have enough timepoints
  if (length(which(!is.na(regiondata))) >= 8) {
    # This will hold a row of features for the gene_region across all groups
    featureVector = c()
    features = c() # A list of labels for the features
    #cat("Processing region",rownames(alldata)[r],"\n")
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
    # Plot to file
    png(paste(imgdir,"/",label,".png",sep=""),width=600,height=300,units = "px", pointsize = 12)
    plot(regiondata,type="l",col="orange",lwd=5,main=label,xaxt="n",ylab="rna-seq expression")
    axis(1,labels=groups$all,at=seq(1,length(groups$all)))
    points(interpolated, col = 2, pch = "*")
    dev.off()
    # Here are features that would 100% be redundant to do for each group
    group = names(groups[1])
    timepoints = which(groups[[group]] %in% groups$all)
    ty = interpolated$y[timepoints]
    tx = interpolated$x[timepoints]
    # Number of peaks
    peaks = length(tx[which(diff(sign(diff(ty)))==-2)])
    featureVector = c(featureVector,peaks)
    features = c(features,paste("number_peaks",group,sep="_"))
    # If each timepoint has a peak
    peakvector = rep(0,length(timepoints))
    names(peakvector) = groups$all
    peakvector[tx[which(diff(sign(diff(ty)))==-2)]] = 1
    featureVector = c(featureVector,peakvector)
    features = c(features,paste(gsub(" ","",names(peakvector)),"haspeak",group,sep="_"))
    # Now iterate through groupings to extract features from interpolated
    for (g in 1:length(groups)){
      #cat("Processing group",g,"of",length(groups),"\n")
      # First split timeseries into groups based on the ages
      group = names(groups[g])
      timepoints = which(groups[[group]] %in% groups$all)
      ty = interpolated$y[timepoints]
      tx = interpolated$x[timepoints]
      # temporal_ac	 The autocorrelation of the timeseries at lag of one for different timepoints
      temporal_ac = acf(ty,plot=FALSE)
      temporal_ac = as.numeric(temporal_ac$acf)[2] 
      featureVector = c(featureVector,temporal_ac)
      features = c(features,paste("temporal_ac_lag1",group,sep=""))
      # temporal_bins   The power of the timeseries collected in different bins to fit the range of expression values
      breaks = c(0,0.1,0.5,1.0,5,10,25,50,75,100,500,1000,5000,10000,20000,30000,40000,50000,300000)
      histy = hist(ty,breaks=breaks,plot=FALSE)
      featureVector = c(featureVector,histy$counts)
      bins = c("histogram_0_0.1","histogram_0.1_0.5","histogram_0.5_1","histogram_1_5","histogram_5_10",
               "histogram_10_25","histogram_25_50","histogram_50_75","histogram_75_100","histogram_100_500",
               "histogram_500_1000","histogram_1000_5000","histogram_5000_10000","histogram_10000_20000",
               "histogram_20000_30000","histogram_30000_40000","histogram_40000_50000","histogram_50000_300000")
      bins = paste(bins,group,sep="_")
      features = c(features,bins)
      # min, max, mean, and sd
      minmaxmean = c(min(ty),max(ty),mean(ty),sd(ty))
      featureVector = c(featureVector,minmaxmean)
      features = c(features,paste(c("min","max","mean","sd"),group,sep="_"))
      # Difference between max and min, min and mean, and mean and max
      tmp = c(abs(max(ty)-min(ty)),abs(mean(ty)-min(ty)),abs(max(ty)-mean(ty)))
      featureVector = c(featureVector,tmp)
      features = c(features,paste(c("abs_diff_max_min","abs_diff_mean_min","abs_diff_mean_max"),group,sep="_"))
      # Entropy of the timeseries
      shannon = shannon.entropy(ty)
      featureVector = c(featureVector,shannon)
      features = c(features,paste("shannon_entropy",group,sep="_"))
      # Kurtosis
      kurt = kurtosis(ty)
      featureVector = c(featureVector,kurt)
      features = c(features,paste("kurtosis",group,sep="_"))
      # Change between timepoints (mean, min, and max)
      change = c(mean(abs(diff(ty))),min(abs(diff(ty))),max(abs(diff(ty))))
      featureVector = c(featureVector,change)
      features = c(features,paste(c("mean_change","min_change","max_change"),group,sep="_"))
      # Skewness
      skew = skewness(ty)
      featureVector = c(featureVector,skew)
      features = c(features,paste("skewness",group,sep="_"))           
    }
    cat("Length of featureVector is",length(featureVector),"\n")
    temporal = rbind(temporal,featureVector)
    collabels = features
  }
}

rownames(temporal) = rowlabels
colnames(temporal) = collabels
ts = list(features=temporal,columns=collabels,rows=rowlabels)
save(ts,file=outfile)