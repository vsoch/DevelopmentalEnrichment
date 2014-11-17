# This script will parse the developmentaEnrichment results, and
# calculate corrected p values (correcting for all tests)

# disorderEnrichment Parsing
setwd("/scratch/users/vsochat/DATA/ALLEN/UberDisorder/developmentRegions1/")
files = list.files(pattern="_developmentRegionEnrichment.Rda")
raw = c()

# Save all Pvalues under .01
numtests = 0
for (f in 1:length(files)){
  cat("Loading file",f,"of",length(files),"\n")
  load(files[f])
  numtests = numtests + nrow(fisher)
  #fisher = fisher[which(fisher$PVALUE <= .05),]
  raw = rbind(raw,fisher)
}

# Let's do FDR correction because we have all the pvalues
FDR = p.adjust(raw$PVALUE,method="fdr")
fisher = as.data.frame(raw,stringsAsFactors=FALSE)
fisher = cbind(fisher,FDR)
fisher$PVALUE= as.numeric(fisher$PVALUE)

# Let's also try bonferroni correction - divide the threshold by the number of tests
rate01 = .01 / numtests
rate05 = .05 / numtests

# This gives us our new rate, how many pass?
bonferroni01 = fisher[which(fisher$PVALUE <= rate01),]
bonferroni05 = fisher[which(fisher$PVALUE <= rate05),]
fdr05 = fisher[which(fisher$FDR <= 0.05),]
fdr01 = fisher[which(fisher$FDR <= 0.01),]

# Here are significant results!
DevelRegions = list(bonferroni01=bonferroni01,bonferroni05=bonferroni05,rate01=rate01,rate05=rate05,fdr05=fdr05,fdr01=fdr01,raw=raw)
save(DevelRegions,file="/scratch/users/vsochat/DATA/ALLEN/UberDisorder/developmentRegions1Results/develRegions_sigresult.Rda")

# How many unique disorders?
length(unique(bonferroni01$DISORDER))
length(unique(bonferroni05$DISORDER))
length(unique(fdr01$DISORDER))
length(unique(fdr05$DISORDER))

cat(unique(bonferroni01$DISORDER),sep="\n",file="/scratch/users/vsochat/DATA/ALLEN/UberDisorder/developmentRegions1Results/uniqueDevelSignificantResult01.txt")
cat(unique(bonferroni05$DISORDER),sep="\n",file="/scratch/users/vsochat/DATA/ALLEN/UberDisorder/developmentRegions1Results/uniqueDevelSignificantResult05.txt")

# For each disorder, we should now make a timeseries of brain maps that show where the enrichment is! (and how it changes)
setwd("/scratch/users/vsochat/DATA/ALLEN/UberDisorder/developmentRegions1Results")
library("scatterplot3d")

# For each disorder
uniquedisorders = as.character(unique(bonferroni01$DISORDER))
uniqueregions = as.character(unique(raw$REGION))

# Do disorders have more than one timepoint?
ridx = c()
for (r in 1:length(uniquedisorders)){
  disorder = uniquedisorders[r]
  tmp = bonferroni01[which(bonferroni01$DISORDER == disorder),]
  timepoints = unique(tmp$TIMEPOINT)
  if (length(timepoints)>1){
    printtime = paste(timepoints,collapse="|")
    cat("Disorder",disorder,"has",length(timepoints),"timepoints:",printtime,"\n")
    ridx = c(ridx,r)
  }
}


# First let's make a matrix of disorder by region
regionmatrix = array(0,dim=c(length(uniquedisorders),length(uniqueregions)))
colnames(regionmatrix) = sort(uniqueregions)
rownames(regionmatrix) = sort(uniquedisorders)

# Parse the result
for (r in 1:nrow(bonferroni01)){
  cat("Processing",r,"of",nrow(bonferroni01),"\n")
  region = as.character(bonferroni01[r,"REGION"])
  disorder = as.character(bonferroni01[r,"DISORDER"])
  regionmatrix[disorder,region] = 1
}

save(regionmatrix,file="Devel_sigresult_bonf01_regionmatrix.Rda")
pdf("regionmatrix_clustering.pdf")
heatmap(regionmatrix)
dev.off()


# Visualize the results
sink(file=paste('DevelSigresult_bonf01_Brainmaps.txt',sep=""))
# Do disorders have more than one timepoint?
for (r in 1:length(uniquedisorders)){
  #cat("Processing",r,"of",nrow(bonferroni01),"\n")
  disorder = gsub(" ","",uniquedisorders[r])
  tmp = bonferroni01[which(bonferroni01$DISORDER == disorder),]
  timepoints = unique(tmp$TIMEPOINT)

  cat("\n",disorder,"\n")
  # Now let's visualize the coordinates
  # If we have more than one timepoint
  if (length(timepoints)>1){
    #par(mfrow = c(1, length(timepoints)))
    for (t in 1:length(timepoints)){
      tp = timepoints[t]
      cat("  ",tp,"\n")
      tmp2 = tmp[which(tmp$TIMEPOINT==tp),]
      re = unique(tmp2$REGION)
      for (rr in 1:length(re)){
        cat("    ",re[rr],"\n")
      }
      #coords = cbind(tmp2$mni_x,tmp2$mni_y,tmp2$mni_z)
      #scd = scatterplot3d(coords,pch=19,color="orange",main=paste("Enriched Regions for",disorder,tp),xlab="MNIx",ylab="MNIy",zlab="MNIz",ylim = c(-100,70),xlim = c(-70,70),zlim = c(-70,80))
      # We need to convert to 2d point labels
      #twodcoords = scd$xyz.convert(coords)
      #text(twodcoords$x,twodcoords$y,labels=tmp2$REGION,cex=.5, pos=4)  
    }
   # Otherwise just do one plot
   } else {
    tmp2 = tmp[which(tmp$TIMEPOINT==timepoints),]
    re = unique(tmp2$REGION)
    cat("  ",timepoints,"\n")
       for (rr in 1:length(re)){
        cat("    ",re[rr],"\n")
   }   
#par(mfrow = c(1,1))
    # Now let's visualize the coordinates
    #coords = cbind(tmp$mni_x,tmp$mni_y,tmp$mni_z)
    #scd = scatterplot3d(coords,pch=19,color="orange",main=paste("Enriched Regions for",disorder,timepoints[1]),xlab="MNIx",ylab="MNIy",zlab="MNIz",ylim = c(-100,70),xlim = c(-70,70),zlim = c(-70,80))
    # We need to convert to 2d point labels
    #twodcoords = scd$xyz.convert(coords)
    #text(twodcoords$x,twodcoords$y,labels=tmp$REGION,cex=.5, pos=4)  
  }
}h
sink()
