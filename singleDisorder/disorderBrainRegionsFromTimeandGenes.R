# Disorder region enrichment
# Given we know points in time when a disorder is manifested
# Given we know genes involved
# Can we pinpoint brain regions likely to be involved with the manifestation of disorder?
# Can we find patterns of expression for genes associated with a disorder?

load("timeseries_featurematrix.Rda")
load("timeseries_genes_list.Rda")
load("genes_47808_timeseries.Rda")
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

# Now let's find these genes in the developmental data
# Which ones are missing?
missing = asd_genes[-which(asd_genes %in% tgenes)]

# Manually fix gene names, yuck
missing[10] = "HOXB"
missing[11] = "HOXD"
missing[18] = "PCDHA"

# Do we want to add these? Maybe not...
contenders = c()
for (g in missing){
  tmp = tgenes[grep(g,tgenes)]
  cat("GENE IS",g,"\nCONTENDERS:",tmp,"\n")
  contenders = c(contenders,tmp)
}

final_genes = interpolated$genes[which(interpolated$genes %in% asd_genes)]

# Before we subset the data, let's interpolate to get missing points
# for each gene.

# Get the genes! Here is the raw timeseries data
subset = interpolated$timeseries[which(interpolated$genes %in% asd_genes),]
ts = list(timeseries=subset,genes=final_genes)
save(ts,file="/scratch/users/vsochat/DATA/ALLEN/NeuroDisorder/genelists/Autism_timeseries_subset_446.Rda")

# Let's take a look at the timecourses
load("/home/vanessa/Documents/Dropbox/Code/R/DevelopmentEnrichment/data/Autism_timeseries_subset_446.Rda")

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
pdf("asd_timeseries_zscore_genes_446.pdf",onefile=TRUE)
for (gene in unique(ts$genes)){
  tmp = Z[wyhich(ts$genes == gene),]
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

# STOPPED HERE!
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
save(asd,file="asd_ts_7136.Rda")

# For each gene, plot across development
pdf("asd_gene_patterns.pdf",onefile=TRUE,width=12)
colors = sample(colours(),length(unique(regions)))
for (gene in names(asd_genes)){
  tmp = asd_genes[[gene]]
  plot(tmp[1,],col=colors[1],pch=19,ylab="rna-seq expression",xaxt="n",main=paste("normalized expression for ",gene),xlab="age")
  axis(1,labels=colnames(tmp),at=seq(1,length(colnames(tmp)))) 
  for (o in 2:nrow(tmp)) {
    points(tmp[o,], col = colors[o],pch=19)
  }
  legend(0,max(tmp),regions,lty=c(1,1),lwd=c(2.5,2.5),col=colors,cex=.6)
}
dev.off()



# For each gene, plot across development

# Question 3: Genes moving together within each timepoint

# Let's try a simple clustering!
heat = pheatmap(Z)


}
dev.off()


stopped here -
get annotations - cell type = something relevant to genes
cluster and label WITH those labels
write about

# Here are the timeseries features
feature

# Iterpolate missing points

# What are the paterns?

names(mala)[grep("infl",tolower(names(mala)))]

