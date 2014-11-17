# Filter down to values / features that we have defined
idx = which(!is.na(gene))
data = gene[idx]

# For each other gene/region in the matrix, compare
# This will be a vector of the number of features we could compare
# We may want to filter this down to only include gene/regions with a reasonable number
ncompared = array(dim=nrow(featurematrix))
names(ncompared) = rownames(featurematrix)

# This will hold the correlation
sim = array(dim=nrow(featurematrix))
names(sim) = rownames(featurematrix)

# This will hold the pvalue
pvalue = array(dim=nrow(featurematrix))
names(pvalue) = rownames(featurematrix)

# This will hold the rsquared
rsquared = array(dim=nrow(featurematrix))
names(rsquared) = rownames(featurematrix)

for (i in 1:length(sim)){
  cat("Processing",i,"of",length(sim),"\n")
  genetocompare = names(sim)[i]
  # This is redundant, but we want to be sure
  ii = which(rownames(featurematrix)==genetocompare)
  rowtocompare = featurematrix[ii,]
  # First filter down to the features we have for our row
  commonfeatures = which(names(rowtocompare) %in% names(data))
  # Now filter down to remove NA
  rowtocompare = rowtocompare[commonfeatures]
  finalcomparison = intersect(which(!is.na(rowtocompare)),which(!is.nan(rowtocompare)))  
  rowtocompare = rowtocompare[finalcomparison]
  # Filter our original down
  datatmp = data[which(names(data) %in% names(rowtocompare))]
  
  if (length(datatmp) > 0) {    
    
    # Save the number of features we are comparing
    ncompared[genetocompare] = length(datatmp)
    
    # Calculate Correlations
    tmp = coef(summary(lm(datatmp~rowtocompare)))
    if ((dim(tmp)[1]>=2) && (dim(tmp)[2]>=1)) { 
      sim[genetocompare]  = coef(summary(lm(datatmp~rowtocompare)))[2,1]
    } else {
      sim[genetocompare] = 0
    }
    # Calculate pvalue and rsquared
    tmp = coef(summary(lm(datatmp~rowtocompare)))
    if ((dim(tmp)[1]>=2) && (dim(tmp)[2]>=4)) { 
      pvalue[genetocompare] = coef(summary(lm(datatmp~rowtocompare)))[2,4]
      rsquared[genetocompare] = summary(lm(datatmp~rowtocompare))$r.squared 
    } else {
      pvalue[genetocompare] = -999
      rsquared[genetocompare] = 0
    }
  } else{
    cat("No common variables for gene",genetocompare,"\n")
    ncompared[genetocompare] = 0
  }
}

# Save result to file     
result = list(rowname=name,pvalue=pvalue,rsquared=rsquared,corr=sim,rowidx=row,numcompared=ncompared)
save(result,file=outfile)