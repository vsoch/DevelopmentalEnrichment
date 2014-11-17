## DEVELOPMENTAL BRAIN
# Find the timepoints / regions in development with enrichment for disorder genes

args <- commandArgs(TRUE)
disorder = args[1]
resultfile = args[2]

load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/rna-genes-52375x525.Rda")
# Here are disorder genes
load("/scratch/users/vsochat/DATA/ALLEN/UberDisorder/Uber5_GenesGr1.Rda")

# This is the time period we will filter to
utero = c("8 pcw","9 pcw","12 pcw","13 pcw","16 pcw","17 pcw","19 pcw","21 pcw","24 pcw","25 pcw","26 pcw","35 pcw","37 pcw")
uidx = rna.genes$columns[which(rna.genes$columns$age %in% utero),]
people = rna.genes$columns[which(rna.genes$columns$age %in% utero),]
probes = rna.genes$rows
probes$gene_symbol =  as.character(probes$gene_symbol)
data = rna.genes$data[,uidx$column_num]

# We are again going to compare to the means across ALL DEVELOPMENT
# Let's calculate the mean for each gene (rows)
genemeans = rowMeans(data)
genesd = apply(data,1,sd)

# We need to immediately eliminate genes with 0 expression - we will add them to 
# the "notweird" counts
upperlimit = genemeans + 3*genesd
lowerlimit = genemeans - 3*genesd
notweirdadditionidx = which(upperlimit == lowerlimit)
genemeans = genemeans[-notweirdadditionidx]
genesd = genesd[-notweirdadditionidx]
upperlimit = upperlimit[-notweirdadditionidx]
lowerlimit = lowerlimit[-notweirdadditionidx]

# Get rid of these genes from the data
data = data[-notweirdadditionidx,]

# Before we filter the probes, save the probes with 0 mean expression (always "not weird")
notweirdaddition = probes$gene_symbol[notweirdadditionidx]
probes = probes[-notweirdadditionidx,]

setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING")
timepoints = utero

# ANALYSIS 3
# This next analysis will look on a regional basis - is the expression in a region 
# of the brain different than other parts of the brain for a particular timepoint?
# This is the same as the first analysis, but we have multiple timepoints
# We would want to see some kind of clustering of regions at a particular timepoint

fisher = c()
cat("Processing disorder",disorder,"\n")
gen = names(genes[[disorder]])
  
for (t in 1:length(timepoints)){
  cat("Processing timepoint",t,"\n")

  # Split into disorder and non disorder genes
  disordergenes = probes$gene_symbol[which(probes$gene_symbol %in% gen)]
  nondisordergenes = probes$gene_symbol[-which(probes$gene_symbol %in% gen)]

  # Get data for that timepoint
  exp = data[,which(people$age %in% timepoints[t])]
  ppl = people[which(people$age %in% timepoints[t]),]

  # Get unique regions for the timepoint
  uniqueregions = as.character(unique(ppl$structure_name))

  # For each region
  for (r in 1:length(uniqueregions)){
    region = uniqueregions[r]

    # Filter data to that
    idx = people$column_num[which(ppl$structure_name == region)]

    if (!is.vector(exp)){
      expregion = exp[,idx]
    } else {
      expregion = exp        
    }
    pplregion = ppl[which(ppl$structure_name == region),]

    # Find genes above or below 3SD of the mean
    weird = probes[c(which(expregion >= upperlimit),which(expregion <= lowerlimit)),]
    notweird = probes[intersect(which(expregion <= upperlimit),which(expregion >= lowerlimit)),]

    # Make a 2x2 table for the regions/disorder
    # High or Low Expression and disorder gene
    a = length(which(unique(weird$gene_symbol) %in% disordergenes))
    # Not high or low expression and disorder gene
    b = length(which(unique(notweird$gene_symbol) %in% disordergenes)) + length(which(unique(notweirdaddition) %in% disordergenes))
    # High or Low Expression and non-disorder gene
    c = length(which(unique(weird$gene_symbol) %in% nondisordergenes))
    # Not high or low expression and non-disorder gene
    d = length(which(unique(notweird$gene_symbol) %in% nondisordergenes)) + length(which(unique(notweirdaddition) %in% nondisordergenes))

    # Fisher's exact test'
    tabley = matrix(c(a,c,b,d),nrow=2,dimnames=list(c("DisorderGenes","~DisorderGenes"),c("WeirdExpression","~WeirdExpression")))

    ft = fisher.test(tabley, conf.level = 0.95)
    fisher = rbind(fisher,c(disorder,timepoints[t],region,ft$p.value))
  }
}
fisher = as.data.frame(fisher,stringsAsFactors=FALSE)
colnames(fisher) = c("DISORDER","TIMEPOINT","REGION","PVALUE")
fisher$PVALUE = as.numeric(fisher$PVALUE)
# Save fisher result to outfile
save(fisher,file=resultfile)
