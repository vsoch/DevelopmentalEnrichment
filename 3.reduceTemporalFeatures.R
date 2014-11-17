# This script will first explore the distributions of timeseries data, and then
# we need to reduce the regional timeseries for a gene to a smaller set

library(corrplot)

setwd("/home/vanessa/Documents/Dropbox/Code/R/DisorderEnrichment/")
source("3.factorAnalysisFunctions.R")

setwd("/home/vanessa/Documents/Work/DEVELOPING")
# Load timeseries
load("timeseries_featurematrix.Rda")
# Get names of genes
load("timeseries_genes_list.Rda")

# First let's look at the distributions for ALL our features!
pdf("timeseries_feature_histograms.pdf",onefile=TRUE)
for (col in 1:ncol(featurematrix)){
  hist(featurematrix[,col],col=sample(colours(),1),main=colnames(featurematrix)[col],xlab="")
}
dev.off()

# Do we have the computational power to look at correlations across all regions/genes?
# Which features are useless? (eg, colSums==0)
tmp = featurematrix[,-which(colSums(featurematrix)==0)] # we lose 4 features - which ones?
colnames(featurematrix)[which(colSums(featurematrix)==0)]

# Find genes with NA (need to troubleshoot)
findna = which(is.na(featurematrix),arr.ind=TRUE)
narows = unique(findna[,1])
nacols = unique(findna[,2])
colnames(featurematrix)[nacols]

tmp[is.na(tmp)] = 0 # Set NA to zero to deal with now
corrall = cor(tmp)
corrplot(corrall, order="hclust",method = "color",tl.pos="n")

# It looks like we can get NA for shannon entropy, temporal_ac, kurtosis, and skewness
# skewness:
# kurtosis:
# temporal_ac:
# shannon_entropy
# Most of these are based on division by 0

# Are any genes missing?
missing = c()
for (g in 1:length(genes)){
  cat("Processing",g,"of",length(genes),"\n")
  gene = genes[g]
  idx = grep(gene,rownames(featurematrix)) 
  if (length(idx)==0){
    missing = c(missing,gene)
  }
}

# Do Factor Anaysis
# NOT DONE YET.  TRYING UNSUPERVISED CLUSTERING FIRST
#pdf("TimeseriesFeatureCorrelations.pdf",onefile=TRUE)
#nagenes = c()
for (g in 1:length(genes)){
  cat(g,"of",length(genes),"\n")
  gene = genes[g]
  idx = grep(gene,rownames(featurematrix)) 
  # Filter down to those rows
  if (length(idx) > 0) {
    tmp = featurematrix[idx,]
    # We need to eliminate columns (features) with all zero
    tmp = tmp[,-which(colSums(tmp)==0)]
   
    ev = eigen(cor(tmp)) # get eigenvalues
    ap = parallel(subject=nrow(tmp),var=ncol(tmp),rep=100,cent=.05)
    nS = nScree(x=ev$values, aparallel=ap$eigen$qevpea)
    plotnScree(nS,main=paste("Scree test for",gene))
    
    # Now do factor analysis
    # Do we want to normalize the weight matrix to Z scores and choose some threshold?
    K = nS$Components$noc
    cat("Number of optimal factors:",K,"\n")
    
    # This won't work if singular
    try {
      fa = factanal(tmp, factors=K, rotation="promax", scores="regression")  
    } except   
    # This will do smoothing
    fa = fa(r=tmp, nfactors=K, rotate="varimax", SMC=FALSE, fm="minres")
    
    # Extracting factor loadings 
    loadings <- loadings(fa)
    # Visualize factor loadings before thresholding
    pheatmap(loadings,main=paste("Loadings before thresholding for",gene))
    # Plot to identify specific variables for each factor
    fa.diagram(loadings,rsize=0.5,cut=0.1,,main=paste("Factor Loadings",gene))
    # Binary matrix to identify the questions pointed in fa.diagram
    aux = likefadiagram(loadings)
    bin.load = aux$binary
    load = aux$load; load.NA = aux$load.NA
    # We don't want to lose our labels
    rownames(aux$load) = attr(aux$load,"dimnames")[[1]]
    colnames(aux$load) = attr(aux$load,"dimnames")[[2]]
    rownames(aux$load.NA) = attr(aux$load.NA,"dimnames")[[1]]
    colnames(aux$load.NA) = attr(aux$load.NA,"dimnames")[[2]]
    rownames(bin.load) = attr(aux$load,"dimnames")[[1]]
    colnames(bin.load) = attr(aux$load,"dimnames")[[2]]
    # Plot new heatmap
    pheatmap(bin.load,cluster_rows=T,cluster_cols=F,main=paste("Final question assignments for",gene))
    
    # Save data, fa model, and factors to file
    readme = c("fa: the factor analysis of data\nF: the matrix F (N people x K factors) to be used as features\nscree: to determine number of factors\nK: number of factors\nquestions: assigned questions to each factor\ndata: raw questionnaire data with NA and missing subjects removed\nloadings: the matrix W (K factors x P questions)\nMatrix was singular, and smoothed for FA")
    tmp = list(fa=fa,F=fa$scores,questionnaire=label,scree=nS,K=K,questions=bin.load,data=data,loadings=loadings,README=readme)
    label
    factorAnalysis[[label]] = tmp
    save(factorAnalysis,file="FactorAnalysisBattery.Rda")
    
  }    
}

# How to do correlation
#corr = cor(tmp)
#corrplot(corr, order="hclust",method = "color",tl.pos="n",sub=paste("Correlation Matrix Gene-region Features for",gene))