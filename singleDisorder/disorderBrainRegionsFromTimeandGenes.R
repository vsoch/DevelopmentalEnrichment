# Disorder region enrichment
# Given we know points in time when a disorder is manifested
# Given we know genes involved
# Can we pinpoint brain regions likely to be involved with the manifestation of disorder?
# Can we find patterns of expression for genes associated with a disorder?

load("timeseries_featurematrix.Rda")
load("timeseries_genes_list.Rda")
load("genes_47808_timeseries.Rda")
load("interpolated_timeseries.Rda")
tgenes = genes

load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/interpolated_timeseries.Rda")
# Let's load the disorder gene lists
load("/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/genelists/malaDisorderGenes.Rda")

# Find autism
contenders = grep("aut",tolower(names(mala)))
names(mala)[contenders]

# Oh, I see autism is second in the list!
disorder = names(mala)[contenders[2]]
gen = mala[[disorder]]

asd_genes = c()
for (g in gen) {
  asd_genes = c(asd_genes,strsplit(g,"[(]")[[1]][1])
}


# Here are highly validated ASD genes from 
# http://www.nature.com/nature/journal/v515/n7526/fig_tab/nature13772_T1.html
# (look in supplementary tables for 108 genes)
# http://www.nature.com/nature/journal/v515/n7526/full/nature13908.html#tables
asd_genes = c("ADNP","ANK2","ARID1B","CHD8","CUL3","DYRK1A","GRIN2B","KATNAL2","POGZ","SCN2A","SUV420H1","SYNGAP1","TBR1","ASXL3","BCL11A","CACNA2D3","MLL3","ASH1L","CTTNBP2","GABRB3","PTEN","RELN","APH1A","CD42BPB","ETFB","NAA15","MYO9B","MYT1L","NR3C2","SETD5","TRIO","MIB1","VIL1","SCN2A","SYNGAP1","CHD8","ARID1B","ANK2","SUV420H1","DYRK1A","GRIN2B","ADNP","TBR1","POGZ","CUL3","KATNAL2","BCL11A","CACNA2D3","MIB1","GABRB3","KMT2C","PTEN","RELN","ASXL3","MYO9B","ETFB","ASH1L","TRIO","NAA15","MYT1L","NR3C2","APH1A","VIL1","CDC42BPB","BIRC6","GALNTL4","MTMR12","EP400","ZNF774","DPP3","FCRL6","SMURF1","TTLL3","UTP6","CARKD","DNAH10","PPM1D","JADE2","S100G","GSE1","CSNK1E","KIAA0100","CACNA1D","RAB2A","PCOLCE","MYOC","SLCO1B1","SLC6A1","ATP1B1","SCARA3","KRT34","BRWD1","PRPF39","SIX2","CCSER1","GGNBP2","NRXN1","WHSC1","KLC1","RANBP17","LEO1","PTPRM","KDM4B","SETBP1","SRPK2","ZNF559","CSDE1","JUP","QRICH1","GSDMC","AGAP2","PLA1A","HDLBP","TGM1","LRRC14","KDM3A","C11orf30","TAF4","TCTE3","CERS4","TCF3","SLCO1B3","CD163L1","NCKAP1","CSTF2T","BRSK2","MYH10","STXBP5","SHANK3","AXL","IQGAP2","UBN2","KDM6B","CAPN12")
# Do we want to add these? Maybe not...
#contenders = c()
#for (g in missing){
#  tmp = tgenes[grep(g,tgenes)]
#  cat("GENE IS",g,"\nCONTENDERS:",tmp,"\n")
#  contenders = c(contenders,tmp)
#}

final_genes = interpolated$genes[which(interpolated$genes %in% asd_genes)]

# Before we subset the data, let's interpolate to get missing points
# for each gene.

# Get the genes! Here is the raw timeseries data
subset = interpolated$timeseries[which(interpolated$genes %in% asd_genes),]
ts = list(timeseries=subset,genes=final_genes)
save(ts,file="/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/genelists/Autism_timeseries_subset_446.Rda")

# Let's take a look at the timecourses
load("/home/vanessa/Documents/Dropbox/Code/R/DevelopmentEnrichment/data/Autism_timeseries_subset_100.Rda")

# Let's get annotations for these genes - it would be helpful to break the
# into groups first
functional_annot = read.csv("/home/vanessa/Documents/Dropbox/Code/R/DevelopmentEnrichment/data/asd100/functional_annotation.tsv",sep="\t",head=TRUE)

# KEGG
functional_annot$KEGG_PATHWAY
functional_annot$OMIM_DISEASE
functinoal_annot$SP_PIR_KEYWORDS
# Think about how to include this

# Libraries for visualizations, scree plot
library(qgraph)
library(psych)
library(pheatmap)
library(nFactors)
library(Hmisc)

# Since we want to look at patterns of expression first (and we will assess levels after) - let's convert the rows into Z scores
convertToZ = function(row) {
  zscore =  (row - mean(row))/sd(row)
  return(zscore)
}
Z = apply(ts$timeseries,2,convertToZ)

# Come up with 16 unique colors for the brain regions
colors = sample(colors(),16)

# Ok, first let's plot the normalized timeseries for each gene
setwd("/home/vanessa/Documents/Dropbox/Code/R/DevelopmentEnrichment/data")
pdf("asd_timeseries_noz_genes_100.pdf",onefile=TRUE)
for (gene in unique(ts$genes)){
  tmp = Z[which(ts$genes == gene),]
  plot(tmp[1,],type="l",col=colors[1],lwd=5,ylab="rna-seq expression",xaxt="n",main=paste("normalized expression for ",gene),xlab="age")
  axis(1,labels=colnames(Z),at=seq(1,length(colnames(Z))))
  for (o in 2:nrow(tmp)){
    lines(tmp[o,], col = colors[o])
  }
  # Now add a legend
  regions = gsub(paste(gene,"_",sep=""),"",rownames(tmp))
  legend(0,max(tmp),regions,lty=c(1,1),lwd=c(2.5,2.5),col=colors,cex=.6)
}
dev.off()


# Questions from the above:
# What (different) genes move together within a single brain region?
# What single genes are the same across all brain regions (and why??)
# What genes move together within a single timepoint during development?
# I need methods to answer each of these questions.

# Question 1: Genes moving together within single brain regions
# First, create separate matrices for each region with Z scores
# What are unique brain regions?
asd_regions = list() 
for (a in 1:length(regions)){
  asd_regions[regions[a]] = c()
}
for (zz in 1:nrow(Z)){
  cat("Processing",zz,"of",nrow(Z),"\n")
  tmp = Z[zz,]
  gene = ts$genes[zz]
  regionname = gsub(paste(gene,"_",sep=""),"",rownames(Z)[zz])
  asd_regions[[regionname]] =  rbind(asd_regions[[regionname]],tmp)
  rownames(asd_regions[[regionname]])[dim(asd_regions[[regionname]])[1]] = gene
}

# Save our progress
asd = list(timeseries_raw = ts$timeseries,genes=ts$genes,timeseries_z = Z,regions_z = asd_regions)
save(asd,file="asd_ts_7136.Rda")

# Now for each region, let's find genes moving together - first plot
pdf("asd_regional_patterns.pdf",onefile=TRUE,width=12)
# An arbitrary value for labeling gene names in plot on xaxis
xlabs = seq(0,31,by=0.07)
colors = sample(colours(),length(unique(ts$genes)))
for (region in names(asd_regions)){
  tmp = asd_regions[[region]]
  plot(tmp[1,],type="l",col=colors[1],lwd=5,ylab="rna-seq expression",xaxt="n",main=paste("normalized expression for ",region),xlab="age",ylim=c(min(tmp),1))
  axis(1,labels=colnames(tmp),at=seq(1,length(colnames(tmp)))) 
  for (o in 2:nrow(tmp)) {
    lines(tmp[o,], col = colors[o])
    text(xlabs[o],tmp[o,1],rownames(tmp)[o],col=colors[o])
  }
}
dev.off()

# WHOA.  Ok, let's try to find groups within each region:
# Now for each region, let's find genes moving together - first heatmaps
pdf("asd_regional_heatmaps.pdf",onefile=TRUE,width=12)
# An arbitrary value for labeling gene names in plot on xaxis
for (region in names(asd_regions)){
  tmp = asd_regions[[region]]
  heatmap(tmp,main=paste("Which timepoints have similar expression for region",region,"in asd?"))
}
dev.off()

# Now let's try MDS
pdf("asd_regional_mds.pdf",onefile=TRUE)
# An arbitrary value for labeling gene names in plot on xaxis
for (region in names(asd_regions)){
  tmp = asd_regions[[region]]
  d = dist(tmp) # euclidean distances between the rows
  fit = cmdscale(d,eig=TRUE, k=2) # k is the number of dim
  fit # view results
  # plot solution
  x = fit$points[,1]
  y = fit$points[,2]
  #z = fit$points[,3]
  plot(x,y,main=paste("Metric MDS for Region",region),pch=19,col=colors)
  #scatterplot3d(x, y,z, xlab="Coordinate 1", ylab="Coordinate 2",zlab="Coordinate4 3",main=paste("Metric MDS for Region",region),col=colors,pch=19)
  text(x, y, labels = row.names(tmp), cex=0.7) 
}
dev.off()

# Interesting - there are genes (APOE, SNAP25) that are different from the rest within a single region across development.  This MDS is with euclidean distance, and we could arguably use something more specific to temporal data.

pdf("asd_regional_clustering.pdf",onefile=TRUE,width=12)
# An arbitrary value for labeling gene names in plot on xaxis
for (region in names(asd_regions)){
  tmp = asd_regions[[region]]
  # Let's find optimal number of clusters
  intern = clValid(tmp, 5:10, clMethods=c("hierarchical","kmeans","pam"),validation="internal")
  optimal = intern@neighbSize     
  # Let's do clustering, cut the tree, and look at ts:
  hc = hclust(dist(tmp))
  groups = cutree(hc,h=optimal)
  for (g in unique(groups)){
    groupmembers = tmp[which(groups==g),]
    groupcolors = colors[which(groups==g)]
    
    # If we have more than one member in the group
    if (!is.vector(groupmembers)) {
      plot(groupmembers[1,],type="l",col=groupcolors[1],lwd=5,ylab="rna-seq expression",xaxt="n",main=paste("normalized expression for ",region,"group",g),xlab="age",ylim=c(min(groupmembers),max(groupmembers)))
     text(xlabs[1],groupmembers[1],rownames(groupmembers)[1],col=groupcolors[1])
     for (o in 2:nrow(groupmembers)) {
      lines(groupmembers[o,], col = groupcolors[o])
      text(xlabs[o],groupmembers[o,1],rownames(groupmembers)[o],col=groupcolors[o])
     }
    axis(1,labels=colnames(groupmembers),at=seq(1,length(colnames(groupmembers)))) 
    
    # If we only have one guy!
    } else {
     plot(groupmembers,type="l",col=groupcolors[1],lwd=5,ylab="rna-seq expression",xaxt="n",main=paste("normalized expression for ",region,"group",g,"gene",rownames(tmp)[which(groups==g)]),xlab="age",ylim=c(min(groupmembers),max(groupmembers)))
    text(xlabs[1],rownames(tmp)[which(groups==g)]
,col=groupcolors[1])
    axis(1,labels=names(groupmembers),at=seq(1,length(names(groupmembers)))) 
   }  
 }
}
dev.off()


# Question 2: Genes that are consistent across all brain regions
# Goal: to make a plot across time that shows gene expression for each gene at timepoint
# Let's break apart data by gene!
uniquegenes = unique(ts$genes)
asd_genes = list()
for (u in uniquegenes) {
  asd_genes[[u]] = Z[which(ts$genes==u),]
}

# Save our progress
asd = list(timeseries_raw = ts$timeseries,genes=ts$genes,timeseries_z = Z,regions_z = asd_regions, genes_z=asd_genes, regions=regions)
save(asd,file="asd_ts_1600.Rda")

# For each gene, plot across development
pdf("asd_gene_patterns_noz.pdf",onefile=TRUE,width=12)
colors = sample(colours(),length(unique(regions)))
for (gene in names(asd_genes)){
  tmp = asd_genes[[gene]]
  plot(tmp[1,],col=colors[1],pch=19,ylab="rna-seq expression",xaxt="n",main=paste("raw expression for ",gene),xlab="age")
  axis(1,labels=colnames(tmp),at=seq(1,length(colnames(tmp)))) 
  for (o in 2:nrow(tmp)) {
    points(tmp[o,], col = colors[o],pch=19)
  }
  legend(0,max(tmp),regions,lty=c(1,1),lwd=c(2.5,2.5),col=colors,cex=.6)
}
dev.off()


# What genes move together within a single timepoint during development?

# For each gene, plot across development

# Question 3: Genes moving together within each timepoint
# This is pretty simple - now we will cluster genes within timepoints. In this case we would map region on as color

asd_timepoints = list()
uniquetimepoints = colnames(ts$timeseries)
for (u in uniquetimepoints){
  # Let's make a matrix of regions by genes
  mat = array(dim=c(length(unique(ts$genes)),length(regions)))
  rownames(mat) = unique(ts$genes)
  colnames(mat) = regions
  # Now fill in
  for (t in 1:nrow(mat)){
    gene = rownames(mat)[t]
    mat[gene,] = Z[which(ts$genes==gene),u]
  }
  asd_timepoints[[u]] = mat
}

# Save our progress
asd = list(timeseries_raw = ts$timeseries,genes=ts$genes,timeseries_z = Z,regions_z = asd_regions, genes_z=asd_genes, regions=regions, timepoints = asd_timepoints)
save(asd,file="asd_ts_1600.Rda")

# Now again, let's try different clustering
# Are there regions with similar expression in a timepoint, and how
# does this change over time?

pdf("asd_timepoint_clusters.pdf",onefile=TRUE,width=6)
for (u in uniquetimepoints){
  tmp = asd_timepoints[[u]]
  plot(hclust(dist(t(tmp))),main=paste("Which regions have similar rna-seq for timepoint",u,"?"),xlab="",sub="")
}
dev.off()


# DO we want to try WCGNA?
powers = c(seq(4,10,by=1),seq(12,20,by=2))
powerTables = vector(mode="list",length=1)
powerTables[[1]] = list(data=pickSoftThreshold(t(mat),powerVector=powers,verbose=2)[[2]])

colors = c("black","red")
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit","Mean Connectivity","Median Connectivity","Max Connectivity")

ylim = matrix(NA,nrow=2,ncol=4)
for (col in 1:length(plotCols)){
  ylim[1,col] = min(ylim[1,col],powerTables[[1]]$data[,plotCols[col]],na.rm=TRUE)
  ylim[2,col] = min(ylim[2,col],powerTables[[1]]$data[,plotCols[col]],na.rm=TRUE)
}
