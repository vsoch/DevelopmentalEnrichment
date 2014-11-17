# This script will calculate similarity between regions within a gene, and output a reduced 
# number of rows.  We will combine gene-regions that are similar, and keep separate
# those that are not

library(mclust)

args <- commandArgs(TRUE)
gene = args[1]
outfile = args[2]  # Output file

# Load the feature matrix
load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/timeseries_featurematrix.Rda")
load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/timeseries_genes_list.Rda")

cat("Gene:",gene,"\n")
cat("Outfile:",outfile,"\n")

# Filter matrix to regions for that gene
searchstring = paste("^",gene,"_",sep="")
idx = grep(searchstring,rownames(featurematrix))
data = featurematrix[idx,]

# We need to save which columns are zeros (this is different from NA)
zerocol = which(colSums(data) == 0)
data = data[,-zerocol]

# We may want to filter this down to only include gene/regions with a reasonable number
ncompared = array(dim=c(nrow(data),nrow(data)))
rownames(ncompared) = rownames(data)
colnames(ncompared) = rownames(data)

# This will hold the correlation
sim = array(dim=c(nrow(data),nrow(data)))
rownames(sim) = rownames(data)
colnames(sim) = rownames(data)

# This will hold the pvalue
pvalue = array(dim=c(nrow(data),nrow(data)))
rownames(pvalue) = rownames(data)
colnames(pvalue) = rownames(data)

# This will hold the rsquared
rsquared = array(dim=c(nrow(data),nrow(data)))
rownames(rsquared) = rownames(data)
colnames(rsquared) = rownames(data)

# For each row
for (r in 1:nrow(data)){
  cat("Processing row",r,"of",nrow(data),"\n")
  rowname = rownames(data)[r]
  rowdata = data[r,]
  
  # Filter down to values / features that we have defined
  idx = intersect(which(!is.na(rowdata)),which(!is.nan(rowdata)))
  rowdata = rowdata[idx]
  
  # Now compare to other regions for the gene
  for (i in 1:nrow(data)){
   cat("Processing",i,"of",nrow(data),"\n")
   regiontocompare = rownames(data)[i]
   if (i != r) {
      # This is redundant, but we want to be sure
      ii = which(rownames(data)==regiontocompare)
      rowtocompare = data[ii,]
      # First filter down to the features we have for our row
      commonfeatures = which(names(rowtocompare) %in% names(rowdata))
      # Now filter down to remove NA
      rowtocompare = rowtocompare[commonfeatures]
      finalcomparison = intersect(which(!is.na(rowtocompare)),which(!is.nan(rowtocompare)))  
      rowtocompare = rowtocompare[finalcomparison]
      # Filter our original down
      datatmp = rowdata[which(names(rowdata) %in% names(rowtocompare))]

      if (length(datatmp) > 0) {      
        # Save the number of features we are comparing
        ncompared[rowname,regiontocompare] = length(datatmp)
    
        # Calculate Correlations
        tmp = coef(summary(lm(datatmp~rowtocompare)))
        if ((dim(tmp)[1]>=2) && (dim(tmp)[2]>=1)) { 
          sim[rowname,regiontocompare] = coef(summary(lm(datatmp~rowtocompare)))[2,1]
        } else {
          sim[rowname,regiontocompare] = 0
        }
        # Calculate pvalue and rsquared
        tmp = coef(summary(lm(datatmp~rowtocompare)))
        if ((dim(tmp)[1]>=2) && (dim(tmp)[2]>=4)) { 
          pvalue[rowname,regiontocompare] = coef(summary(lm(datatmp~rowtocompare)))[2,4]
          rsquared[rowname,regiontocompare] = summary(lm(datatmp~rowtocompare))$r.squared 
        } else {
          pvalue[rowname,regiontocompare] = 1
          rsquared[rowname,regiontocompare] = 0
        }
      } else {
        cat("No common variables for gene",genetocompare,"\n")
        ncompared[rowname,regiontocompare] = 0
      }
    } else {
        pvalue[rowname,regiontocompare] = 1
        rsquared[rowname,regiontocompare] = 0 
        sim[rowname,regiontocompare] = 0
        ncompared[rowname,regiontocompare] = 0
    }
  }
}

# correct pvalues!
adjusted = array(dim=c(nrow(data),nrow(data)))
adjusted[,] = p.adjust(pvalue,method="fdr")
final = array(0,dim=c(nrow(data),nrow(data)))
final[adjusted <= .01] = 1
diag(final) = 1
rownames(final) = rownames(data)
colnames(final) = rownames(data)

# If we have a model with different groups
if (!all(colSums(final)==ncol(final))) {
  fit = Mclust(final,warn=FALSE)
  groups = fit$classification
  names(groups) = rownames(final)
  
  # Prepare new data
  reduced = array(dim=c(length(unique(groups)),ncol(featurematrix)))
  colnames(reduced) = colnames(featurematrix)

  # If we have unique groups
  if (length(unique(groups))>1) {
    # We will create a list of labels, the regions included in each merged group
    labels = list()
    nregions = c()
    for (u in unique(groups)){
      members = names(groups[groups == u])
      mergedname = paste(gene,"_",u,sep="")
      labels[[mergedname]] = members
      if (length(members)>1) {
        nregions = c(nregions,length(members))
        features = colMeans(data[which(rownames(data)%in%members),])
      } else {
        nregions = c(nregions,1)
        features = data[which(rownames(data) == members),]
      }
      reduced[u,colnames(data)] = features
    }

    # Add our row names
    rownames(reduced) = names(labels)

    # Add zeros to columns that were originally all zero
    reduced[,as.numeric(zerocol)] = 0

    # Calculate weights for each row based on the number of regions
    weights = nregions / nrow(data)
    result = list(gene=gene,reduced=reduced,regions=labels,weights=weights,pvalue=pvalue,rsquared=rsquared,corr=sim,adjusted=adjusted,final=final,clustering=fit,numcompared=ncompared,featurematrixidx = idx,zerocol=zerocol)
    # If we really just have one group
    } else {      
      # Prepare new data
      reduced = array(dim=ncol(featurematrix))
      names(reduced) = colnames(featurematrix)
      
      # We will create a list of labels, the regions included in each merged group
      labels = list()
      members = rownames(data)
      nregions = length(members)
      mergedname = paste(gene,"_",1,sep="")
      labels[[mergedname]] = members
      features = colMeans(data[which(rownames(data)%in%members),])
      reduced[colnames(data)] = features
      # Add zeros to columns that were originally all zero
      reduced[as.numeric(zerocol)] = 0
      # Weights will be 1
      weights = nregions / nrow(data)  
      result = list(gene=gene,reduced=reduced,regions=labels,weights=weights,pvalue=pvalue,rsquared=rsquared,corr=sim,adjusted=adjusted,final=final,clustering=NA,numcompared=ncompared,featurematrixidx = idx,zerocol=zerocol)
    }
  } else {
  # Prepare new data
  reduced = array(dim=ncol(featurematrix))
  names(reduced) = colnames(featurematrix)
  # We will create a list of labels, the regions included in each merged group
  labels = list()
  members = rownames(data)
  nregions = length(members)
  mergedname = paste(gene,"_",1,sep="")
  labels[[mergedname]] = members
  features = colMeans(data[which(rownames(data)%in%members),])
  reduced[colnames(data)] = features
  # Add zeros to columns that were originally all zero
  reduced[as.numeric(zerocol)] = 0
  # Weights will be 1
  weights = nregions / nrow(data)  
  result = list(gene=gene,reduced=reduced,regions=labels,weights=weights,pvalue=pvalue,rsquared=rsquared,corr=sim,adjusted=adjusted,final=final,clustering=NA,numcompared=ncompared,featurematrixidx = idx,zerocol=zerocol)
}

save(result,file=outfile)
