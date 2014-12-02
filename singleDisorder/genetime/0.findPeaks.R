# FindPeaks: extract features to predict time period of disorders

# I give you a gene
# you give me back brain regions and times when it's doing something
# for each result: let me adjust threshold of somerthing, so i can see what happens to other things when i change gene?

# show database of research findings
# show correlation with topic maps

# IDEA
# Can we - create a distribution of times based on set of genes?
# Then assess the distribution as probability, and determine time(s) when
# disorder likely to manifest or change?

# Assumptions: we can look at expression / rna / dna / other data, and determine "important"
# We can understand genes / rna / etc. based on categories of importance - eg, important for function X, in process Y, development of Z.
# If we have probability functions for RELATED disorder genes based on networks, we can have different probability distributios for different behaviors.

# This first script will be a test, very simple - can we extract from the rna data:
# the peak value, for each brain region, 

# How do we know what a high expression value is? Based on region or gene?
# Different ideas to try:

# Absolute max of a region
# Max values in top .05 of distribution
# Peaks anywhere

# TODO: Get other data from Sherlock, for whole story

# For finding peaks
library(quantmod)

setwd("/home/vanessa/Documents/Work/DEVELOPING")
load("timeseries_featurematrix.Rda")
load("timeseries_genes_list.Rda")
load("genes_47808_timeseries.Rda")
load("interpolated_timeseries.Rda")
tgenes = genes

# Come up with 16 unique colors for the brain regions
colors = sample(colors(),16)
colors = sample(rainbow(),16)
line_types = sample(seq(1,6),16,replace=TRUE)

setwd("/home/vanessa/Documents/Dropbox/Code/R/DevelopmentEnrichment/singleDisorder/genetime")


# Come up with 16 unique colors for the brain regions
colors = sample(colors(),16)
colors = rainbow(16)
line_types = sample(seq(1,6),16,replace=TRUE)

# Helper functions
# Set threshold to get top 1%
# All values should be positive, so we don't need to take abs.value
thresh = .95

# Write a function to return threshold for each row
get_top_quantile = function(row) {
qpos = quantile(row,thresh)
  return(qpos)
}

# Interate through ASD 100 genes
asd_genes = c("ADNP","ANK2","ARID1B","CHD8","CUL3","DYRK1A","GRIN2B","KATNAL2","POGZ","SCN2A","SUV420H1","SYNGAP1","TBR1","ASXL3","BCL11A","CACNA2D3","MLL3","ASH1L","CTTNBP2","GABRB3","PTEN","RELN","APH1A","CD42BPB","ETFB","NAA15","MYO9B","MYT1L","NR3C2","SETD5","TRIO","MIB1","VIL1","SCN2A","SYNGAP1","CHD8","ARID1B","ANK2","SUV420H1","DYRK1A","GRIN2B","ADNP","TBR1","POGZ","CUL3","KATNAL2","BCL11A","CACNA2D3","MIB1","GABRB3","KMT2C","PTEN","RELN","ASXL3","MYO9B","ETFB","ASH1L","TRIO","NAA15","MYT1L","NR3C2","APH1A","VIL1","CDC42BPB","BIRC6","GALNTL4","MTMR12","EP400","ZNF774","DPP3","FCRL6","SMURF1","TTLL3","UTP6","CARKD","DNAH10","PPM1D","JADE2","S100G","GSE1","CSNK1E","KIAA0100","CACNA1D","RAB2A","PCOLCE","MYOC","SLCO1B1","SLC6A1","ATP1B1","SCARA3","KRT34","BRWD1","PRPF39","SIX2","CCSER1","GGNBP2","NRXN1","WHSC1","KLC1","RANBP17","LEO1","PTPRM","KDM4B","SETBP1","SRPK2","ZNF559","CSDE1","JUP","QRICH1","GSDMC","AGAP2","PLA1A","HDLBP","TGM1","LRRC14","KDM3A","C11orf30","TAF4","TCTE3","CERS4","TCF3","SLCO1B3","CD163L1","NCKAP1","CSTF2T","BRSK2","MYH10","STXBP5","SHANK3","AXL","IQGAP2","UBN2","KDM6B","CAPN12")

# Filter down to ones we have data for
asd_genes = asd_genes[which(asd_genes %in% interpolated$genes)]

# First - let's make plots of each gene
pdf("asd_bygene_95.pdf",onefile=TRUE,width=12)

for (gene in asd_genes) {

  # Get the subset of the data
  subset = interpolated$timeseries[which(interpolated$genes==gene),]
  regions = gsub(paste(gene,"_",sep=""),"",rownames(subset))

  # First plot the gene
  plot(subset[1,],col=colors[1],pch=19,ylab="rna-seq expression",xaxt="n",main=paste("rna-seq expression for ",gene),xlab="age",ylim=c(min(subset),max(subset)))
  axis(1,labels=colnames(subset),at=seq(1,length(colnames(subset)))) 
  for (o in 2:nrow(subset)) {
    points(subset[o,], col = colors[o],pch=19) 
  }
  legend(0,max(subset),regions,col=colors,cex=.6)

  # Method 1: Maximum value, period
  # Get maximum values for each region (length) labeled by timepoint
  maxvalues = list()
  maxes = apply(subset,1,max)
  for (m in 1:length(maxes)) {
    region = regions[m]
    idx = which(subset[m,]==maxes[m])
    value = subset[m,which(subset[m,]==maxes[m])]
    names(value) = names(idx)
    maxvalues = c(maxvalues,value)
  }

  # Plot the max values
  labels = colnames(subset)
  uniquetimes = unique(names(maxvalues))

  # Plot the first to generate the plot
  time = uniquetimes[1]
  x = which(labels == time)
  values =  maxvalues[which(names(maxvalues) == time)]
  x = rep(x,length(values))  
  r = which(names(maxvalues) == time)
  plot(x,as.numeric(values),col=colors[r],pch=19,ylab="rna-seq expression",xaxt="n",main=paste("timepoints of max values of each region for",gene),xlab="age",ylim=c(min(subset),max(subset)),xlim=c(1,ncol(subset)))
  axis(1,labels=colnames(subset),at=seq(1,length(colnames(subset)))) 

  for (m in 2:length(uniquetimes)) {
    # Get the x coordinate
    time = uniquetimes[m]
    x = which(labels == time)
    # The color is just the value of count, corresponding to row (region)
    r = which(names(maxvalues) == time)  
    values =  maxvalues[which(names(maxvalues) == time)]
    x = rep(x,length(values))
    points(x,as.numeric(values), col = colors[r],pch=19) 
  }
  legend(20,max(subset)/2,lty=c(1,1),regions,col=colors,cex=.6)


  # Method 2: Max values in top .01 of distribution

  # Try getting quantiles for top values
  quantiles = apply(subset,1,get_top_quantile)

  # Get values above threshold for each row
  maxvalues = c()
  timepoints = c()
  reg = c()
  for (r in 1:nrow(subset)) {
    region = r
    idx = which(subset[r,]>=quantiles[r])
    maxvalues = c(maxvalues,subset[r,idx])
    timepoints = c(timepoints,colnames(subset)[idx])
    reg = c(reg,rep(region,length(idx)))
  }

  # Now plot!
  uniquetimes = unique(timepoints)

  # Plot the first to generate the plot
  time = uniquetimes[1]
  x = which(labels == time)
  values =  maxvalues[which(timepoints == time)]  
  x = rep(x,length(values))
  r = reg[which(timepoints==time)]
  plot(x,as.numeric(values),col=colors[r],pch=19,ylab="rna-seq expression",xaxt="n",main=paste("timepoints of top",1-thresh,"of values for each region for",gene),xlab="age",ylim=c(min(subset),max(subset)),xlim=c(1,ncol(subset)))
  axis(1,labels=colnames(subset),at=seq(1,length(colnames(subset)))) 

  for (m in 2:length(uniquetimes)) {
    # Get the x coordinate
    time = uniquetimes[m]
    x = which(labels == time)
    # The color is just the value of count, corresponding to row (region)
    values =  maxvalues[which(timepoints == time)]
    x = rep(x,length(values))
    r = reg[which(timepoints==time)]
    points(x,as.numeric(values), col = colors[r],pch=19)
  }
  legend(20,max(subset)/2,lty=c(1,1),regions,col=colors,cex=.6)

}
dev.off()

# SAVING DATA -------------------------------------------------------------

# We need to also look at change

# Let's save a long list of these important timepoints, we want:
# GENE	TIMEPOINT	REGION	VALUE

# We will use method 2, so we can try two thresholds (.99 and .95)
genetimes = c()
thresh = .95

for (gene in asd_genes) {

  # Get the subset of the data
  subset = interpolated$timeseries[which(interpolated$genes==gene),]
  regions = gsub(paste(gene,"_",sep=""),"",rownames(subset))

  # Method 2: Top thresh%age
  # Try getting quantiles for top values
  quantiles = apply(subset,1,get_top_quantile)

  # If any quantiles == 0, then set to impossibly large number
  quantiles[which(quantiles == 0)] = 99999

  # Get values above threshold for each row
  maxvalues = c()
  timepoints = c()
  reg = c()
  for (r in 1:nrow(subset)) {
    region = r
    idx = which(subset[r,]>=quantiles[r])
    maxvalues = c(maxvalues,subset[r,idx])
    timepoints = c(timepoints,colnames(subset)[idx])
    reg = c(reg,rep(region,length(idx)))
  }

  # Now let's make into a data frame
  result = data.frame(gene=rep(gene,length(maxvalues)),timepoint = timepoints,region=regions[reg],value=maxvalues,thresh=rep(thresh,length(maxvalues)))

   genetimes = rbind(genetimes,result)
}

genetimes = list(matrix=genetimes,timepoints=labels,regions=regions,genes=asd_genes)
save(genetimes,file="genetimes_asd100_99.Rda")

# Now - let's make plots by region!

pdf("asd_byregion_99.pdf",onefile=TRUE,width=12)

for (r in 1:length(regions)) {

  region = regions[r]
  subset = genetimes[which(genetimes$region == region),]
  regioncolor = colors[r]

  # First plot all the genes for a single region
  idx = grep(region,rownames(interpolated$timeseries))
  singleregion = interpolated$timeseries[idx,]
  # Not sure if we want to plot anything here...

  # Plot the gene times in the region
  # Here are the x values
  x = match(subset$timepoint,colnames(interpolated$timeseries))
  plot(x,subset$value,col=regioncolor,pch=19,ylab="rna-seq expression",xaxt="n",main=paste("important timepoints and genes for",region),xlab="age",ylim=c(min(subset$value),max(subset$value)))
  axis(1,labels=colnames(interpolated$timeseries),at=seq(1,length(labels)))
  # Now let's add text for genes
  text(x,subset$value,label=subset$gene,cex=0.6,pos=4) 
}
dev.off()
