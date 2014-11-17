# This script will calculate similarity between meta genes in the reduced temporal feature matrix

bargs <- commandArgs(TRUE)
gene = args[1]
outfile = args[2]  # Output file

# Load the feature matrix
load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/reduced_feature_matrix_51641.Rda")

cat("Gene:",gene,"\n")
cat("Outfile:",outfile,"\n")

# Find our meta gene in the matrix
data = temporal$matrix
idx = which(rownames(data) == gene)

# Grab the complete feature vector
featureVector = data[idx,]
# Get rid of NA values, if there are any
if (any(is.na(featureVector))){
  featureVector = featureVector[-which(is.na(featureVector))]
}

# These vectors will hold a similarity score to all other meta genes
euc.sim = array(dim=nrow(data))
cos.sim = array(dim=nrow(data))
max.sim = array(dim=nrow(data))
man.sim = array(dim=nrow(data))
can.sim = array(dim=nrow(data))
bin.sim = array(dim=nrow(data))
mink.sim = array(dim=nrow(data))
names(euc.sim) = rownames(data)
names(cos.sim) = rownames(data)
names(max.sim) = rownames(data)
names(man.sim) = rownames(data)
names(can.sim) = rownames(data)
names(bin.sim) = rownames(data)
names(mink.sim) = rownames(data)


# Here are functions for euclidean distance, cosine, etc.
euc.dist = function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
cos.dist = function(x1, x2) crossprod(x1, x2)/sqrt(crossprod(x1) * crossprod(x2))

# For each gene, filter down to features that we have in common
for (i in 1:nrow(data)){
  cat("Calculating similarity values for row",i,"\n")
  featureVector2 = data[i,]
  # Filter down to those we have in common
  commonFeatures = which(names(featureVector2) %in% names(featureVector))
  featureVector2 = featureVector2[commonFeatures]
  # Now get rid of NA in feature vector 2
  # If we have any features:
  if (any(is.na(featureVector2))){
    featureVector2 = featureVector2[-which(is.na(featureVector2))]
  }
  # Now limit vector 1 to those
  commonFeatures = which(names(featureVector) %in% names(featureVector2))
  featureVector1 = featureVector[commonFeatures]
  euc.sim[i] = euc.dist(featureVector1,featureVector2)
  cos.sim[i] = cos.dist(featureVector1,featureVector2)
  max.sim[i] = dist(rbind(featureVector1,featureVector2),method="maximum")
  man.sim[i] = dist(rbind(featureVector1,featureVector2),method="manhattan")
  can.sim[i] = dist(rbind(featureVector1,featureVector2),method="canberra")
  bin.sim[i] = dist(rbind(featureVector1,featureVector2),method="binary")
  mink.sim[i] = dist(rbind(featureVector1,featureVector2),method="minkowski")
}

# Save to list
sim = list(gene=gene,rowidx = idx,euclidean=euc.sim,cosine=cos.sim,maximum=max.sim,manhattan=man.sim,canberra=can.sim,binary=bin.sim,minkowski=mink.sim)
save(sim,file=outfile)
