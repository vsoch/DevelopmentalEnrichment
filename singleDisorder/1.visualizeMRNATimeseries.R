# Visualize mrna enrichment across development

# We need to understand how to best separate this data into different time periods,
# and it makes most sense to do this based on patterns of enrichment!

load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/rna-genes-52375x525.Rda")
# Here are disorder genes
load("/scratch/users/vsochat/DATA/ALLEN/UberDisorder/Uber5_GenesGr1.Rda")

# For each brain region, look at change over time
uniqueregions = as.character(unique(rna.genes$columns$structure_name))
timepoints = as.character(unique(rna.genes$columns$age))
uniquegenes = as.character(unique(rna.genes$rows$gene_symbol))


# For each gene
allseries = list()
for (g in 24:length(uniquegenes)){
  cat("Processing gene",g,"of",length(uniquegenes),"\n")
  gene = uniquegenes[g]
  # Create a matrix to hold a timeseries for each region
  timeseries = array(dim=c(length(uniqueregions),length(unique(timepoints))))
  rownames(timeseries) = uniqueregions
  colnames(timeseries) = timepoints
  # Get mean in data for that gene
  # For each region, look at change over time
  for (r in 1:length(uniqueregions)){
    region = uniqueregions[r]
    # This is for multiple timepoints
    cols = rna.genes$columns[which(rna.genes$columns$structure_name %in% region),]
    rows = rna.genes$rows[which(rna.genes$rows$gene_symbol == gene),]
    # Put into the timeseries at the right spots
    ages = as.character(unique(cols$age))
    for (a in ages){
      value = rna.genes$data[rows$row_num[which(rows$gene_symbol==gene)],cols$column_num[which(cols$age==a)]]
      if (!is.vector(value)) {
        value = colMeans(value)
      }
      value = mean(as.numeric(value))
      timeseries[region,a] = value
    }
  }
  allseries[[gene]] = timeseries
}
save(allseries,file="/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/genes_47808_timeseries.Rda")

# Now let's plot the timeseries - work on local machine so we can visualize!
load("/home/vanessa/Documents/Work/DEVELOPING/genes_47808_timeseries.Rda")
load("/home/vanessa/Documents/Work/DEVELOPING/rna-genes-52375x525.Rda")
setwd("/home/vanessa/Documents/Work/DEVELOPING")

# For each gene, take a look at timeseries
genes = names(allseries)
pdf("GeneTimeseriesInAllRegions100.pdf")
for (g in 1:100){
  gene = genes[g]
  timeseries = allseries[[gene]]
  # Only add regions that we have at least 8 points for
  tmp = c()
  for (t in 1:nrow(timeseries)){
    if (sum( !is.na( timeseries[t,] ) ) >= 8){
      tmp = rbind(tmp,timeseries[t,])
    }
  }
  # Plot gene
  # Get max and min
  ymax = max(tmp[!is.na(tmp)])
  ymin = min(tmp[!is.na(tmp)])
  plot(tmp[1,],type="l",col="orange",lty=1,lwd=3,main=gene,ylab="rna-seq expression",xlab="age",xaxt="n",ylim=c(ymin,ymax))
  axis(1,labels=colnames(tmp)[c(1,4,9,14,19,24,29)],at=c(0,5,10,15,20,25,30))
  for (t in 1:nrow(tmp)){
    lines(tmp[t,],col=sample(colours(),1))
  }
}
dev.off()

# Next - finish feature extraction scripts for timeseries.