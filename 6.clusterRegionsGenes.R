# This script will read in the reduced temporal feature matrix and do basic clustering to find groups of genes with similar expression patterns.

setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/finalSimMatrices")
matrices = list.files()

# HIERARCHICAL CLUSTERING -----------------------------------------------------------------------------------
# We will start with standard euclidean distance

load("euclidean_scores.Rda")  # matrix

# STOPPED HERE - need to look at timeseries of groups in clustering!
setwd("/scratch/PI/dpwall/ALLEN-BRAIN/DEVELOPING")
load("timeseries_featurematrix.Rda")

# Here are our 40 groups
# load("euclidean_groups40.Rda")

data = temporal$matrix

# First let's look at the mean features of each group
uniquegroups = unique(groups)
groupMeans = c()
groupids = c()
for (g in uniquegroups){
  members = names(which(groups == as.numeric(g)))
  tmp = data[which(rownames(data)%in%members),]
  if (!is.vector(tmp)){
    groupMeans = cbind(groupMeans,colMeans(tmp,na.rm=TRUE))
    groupids = c(groupids,g)
  } else {
    cat("Group",g,"only has one member!\n")
  }
}
rownames(groupMeans) = colnames(data)
colnames(groupMeans) = groupids

pdf("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/img/group40features.pdf",onefile=TRUE)
for (g in 1:ncol(groupMeans)){
  hist(groupMeans[g,],main=paste("Group Means for feature",rownames(groupMeans)[g]))
}

# Group 27 only has one member!
# Group 30 only has one member!
# Group 31 only has one member!
# Group 32 only has one member!
# Group 33 only has one member!
# Group 34 only has one member!
# Group 36 only has one member!
# Group 37 only has one member!
# Group 38 only has one member!


load("genes_47808_timeseries.Rda")
pdf("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/img/group40tsall.pdf",onefile=TRUE)
# Let's visualize group timeseries
for (t in uniquegroups){
  members = names(which(groups == as.numeric(t)))
  cat("Processing group",t,"of",length(uniquegroups),"\n")
  if (length(members) != 1) {
      # This will hold all timeseries for the group
      allts = c()
      # Get the genes/regions for each!
      for (g in members){
        gene = substring(g,1,nchar(g)-1)
        regions = gsub("_"," ",gsub(gene,"",temporal$rowlookup[[g]]))
        # Now get the original timeseries
        ts = allseries[[gsub("_","",gene)]]
        if (!is.null(ts)){ 
         # Match the rownames
          rownames(ts)  = gsub("[)]|[(]|[/]","",rownames(ts))
          # Get the data!
          tmp = ts[regions,]
          # We don't need rownames for now'
          # Restore row names to original so we don't have duplicates'
          #if (is.vector(tmp)){
          #  rownames(allts)[dim(allts)[1]] = paste(gene,regions,sep="")
          #} else {
          #  rownames(tmp) = paste(gene,rownames(tmp),sep="")
          #}
          allts = rbind(allts,tmp)
         }
      }
    # Now plot the group to file
    ymax = max(allts[!is.na(allts)])
    ymin = min(allts[!is.na(allts)])
    # Plot the mean timeseries for the group
    #meants = colMeans(allts,na.rm=TRUE)
    #plot(meants,type="l",col="orange",lty=1,lwd=3,main=paste("Group",t,"with N =",length(members)),ylab="rna-seq expression",xlab="age",xaxt="n",ylim=c(ymin,ymax))
   # Plot first row to initialize plot
    plot(allts[1,],type="l",col="orange",lty=1,lwd=3,main=paste("Group",t,"with N =",length(members)),ylab="rna-seq expression",xlab="age",xaxt="n",ylim=c(ymin,ymax))
  axis(1,labels=colnames(allts)[c(1,4,9,14,19,24,29)],at=c(0,5,10,15,20,25,30))
  for (t in 2:nrow(allts)){
   lines(allts[t,],col=sample(colours(),1))
  }

  } else {
      cat("Group",t,"only has one member!\n")
  }
   
}

dev.off()


# SELF ORGANIZING MAP -----------------------------------------------------------------------------------
# Let's try creating a self organizing map!
load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/reduced_feature_matrix_51641.v2.Rda")
data = temporal$matrix

# This is slightly irresponsible, but let's replace NaN with 0 for now.
data[is.na(data)] = 0
library('kohonen')

# We are choosing ~40 groups based on clustering above
som = som(data, grid = somgrid(6, 7, "hexagonal"))

# Here are the class assignments for each row in matrix
som$unit.classif

labels = rownames(data)

# Create a vector (with indices corresponding to som$grid$pts) with labels of assigned terms
somLabels = c()
geneGroups = c()
for (t in 1:length(som$unit.classif)){
  tmp = labels[which(som$unit.classif==t)]
  geneGroups = c(geneGroups,list(tmp))
  tmp = paste(tmp,collapse="\n")
  somLabels = c(somLabels,tmp)
}

# Here are the corresponding coordinates of the som
som$grid$pts

# Save to file
geneGrid = list(som=som,nodeLabels=somLabels,nodeGroups=geneGroups,alllabels=labels)
save(geneGrid,file="TemporalSOM42.Rda")

# Let's quickly plot!
plot(geneGrid$som$grid$pts,main="Brainspan Gene Expression Similarity",col="#CCCCCC",xlab="Nodes",ylab="Nodes",pch=19,cex=1)
text(geneGrid$som$grid$pts,geneGrid$nodeLabels,cex=.4)

# clvalid?
# cl = clValid(data, 2:6, clMethods=c("hierarchical","kmeans","pam"),validation="internal")

# VISUALIZATION
# What we need to do now is visualize timeseries for each node - what is the pattern?
# Load original timeseries data!
load("../genes_47808_timeseries.Rda")
genes = names(allseries)
pdf("SOMGroupsMeans_42Node.pdf",onefile=TRUE)
uniquegroups = sort(unique(som$unit.classif))
for (t in uniquegroups){
  group = labels[which(som$unit.classif==t)]

  if (length(group) > 0) {
    # This will hold all timeseries for the group
    allts = c()
    # Get the genes/regions for each!
    for (g in group){
      gene = substring(g,1,nchar(g)-1)
      regions = gsub("_"," ",gsub(gene,"",temporal$rowlookup[[g]]))
      # Now get the original timeseries
      ts = allseries[[gsub("_","",gene)]]
      if (!is.null(ts)){ 
       # Match the rownames
        rownames(ts)  = gsub("[)]|[(]|[/]","",rownames(ts))
        # Get the data!
        tmp = ts[regions,]
        # We don't need rownames for now'
        # Restore row names to original so we don't have duplicates'
        #if (is.vector(tmp)){
        #  rownames(allts)[dim(allts)[1]] = paste(gene,regions,sep="")
        #} else {
        #  rownames(tmp) = paste(gene,rownames(tmp),sep="")
        #}
        allts = rbind(allts,tmp)
       }
    }
  # Now plot the group to file
  ymax = max(allts[!is.na(allts)])
  ymin = min(allts[!is.na(allts)])
  # Plot the mean timeseries for the group
  meants = colMeans(allts,na.rm=TRUE)
  plot(meants,type="l",col="orange",lty=1,lwd=3,main=paste("SOM Group",t,"with N =",length(group)),ylab="rna-seq expression",xlab="age",xaxt="n",ylim=c(ymin,ymax))

  # Plot first row to initialize plot
  #plot(allts[1,],type="l",col="orange",lty=1,lwd=3,main=paste("SOM Group",t,"with N =",length(group)),ylab="rna-seq expression",xlab="age",xaxt="n",ylim=c(ymin,ymax))
  #axis(1,labels=colnames(allts)[c(1,4,9,14,19,24,29)],at=c(0,5,10,15,20,25,30))
  #for (t in 2:nrow(allts)){
  # lines(allts[t,],col=sample(colours(),1))
  #}
  }
}
dev.off()

