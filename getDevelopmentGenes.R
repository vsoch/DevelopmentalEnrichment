## DEVELOPMENTAL BRAIN
# This script will, for each disorder, get the genes with significant result

#args <- commandArgs(TRUE)
#disorder = args[1]
#resultfile = args[2]

load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/rna-genes-52375x525.Rda")
# Here are disorder genes
load("/scratch/users/vsochat/DATA/ALLEN/UberDisorder/Uber5_GenesGr1.Rda")
# Here are the significant results, we will look at bonferroni .01
load("/scratch/users/vsochat/DATA/ALLEN/UberDisorder/developmentRegions1Results/develRegions_sigresult.Rda") # DevelRegions$bonferroni01

# Here is the new results file to write (with disorders, regions, genes)
resultfile = "/scratch/users/vsochat/DATA/ALLEN/UberDisorder/developmentRegions1Results/SigDisorderRegionGenes.Rda"

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
# Important - row names are still the old (before elimination of genes)

# Before we filter the probes, save the probes with 0 mean expression (always "not weird")
notweirdaddition = probes$gene_symbol[notweirdadditionidx]
probes = probes[-notweirdadditionidx,]

setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING")

# ANALYSIS 3
# This next analysis will look on a regional basis - is the expression in a region 
# of the brain different than other parts of the brain for a particular timepoint?
# This is the same as the first analysis, but we have multiple timepoints
# We would want to see some kind of clustering of regions at a particular timepoint

devel = DevelRegions$bonferroni01
uniquedisorders = unique(devel$DISORDER)

# We will save all results here
results = c()

for (d in 1:length(uniquedisorders)){
  disorder = uniquedisorders[d]
  cat("Processing disorder",disorder,d,"of",length(uniquedisorders),"\n")
  gen = names(genes[[disorder]])

  # Filter data down to results for that disorder
  subdevel = devel[which(devel$DISORDER==disorder),]

  # Split into disorder and non disorder genes
  disordergenes = probes$gene_symbol[which(probes$gene_symbol %in% gen)]
  nondisordergenes = probes$gene_symbol[-which(probes$gene_symbol %in% gen)]

  # For each result in devel
  for (e in 1:nrow(subdevel)){
    timepoint = subdevel$TIMEPOINT[e]
    region = subdevel$REGION[e]

    #VANESSA - there should only be ONE PVALUE here.
    pval = subdevel$PVALUE
    corrected = subdevel$FDR
    cat("Processing timepoint",t,",region",region,"\n")

    # Get data for that timepoint
    exp = data[,which(people$age %in% timepoint)]
    ppl = people[which(people$age %in% timepoint),]

    # Filter data to the region
    idx = people$column_num[which(ppl$structure_name == region)]

    if (!is.vector(exp)){
      expregion = exp[,idx]
    } else {
      expregion = exp        
    }
    names(expregion) = rownames(exp)
    pplregion = ppl[which(ppl$structure_name == region),]

    # Find genes above or below 3SD of the mean
    # These are genes that are up/down, regardless of the disorder
    upidx = which(expregion >= upperlimit)
    downidx = which(expregion <= lowerlimit)
    # Here are all upregulated genes
    upgenes = probes[upidx,]
    downgenes = probes[downidx,]
    # Now filter the up and down genes to only those in our disorder
    upgenes = upgenes[which(upgenes$gene_symbol %in% disordergenes),]
    downgenes = downgenes[which(downgenes$gene_symbol %in% disordergenes),]
    # Save the expression values of the up and down genes
    # Important! We are indexing based on the row name here, not index
    upgenesexp = expregion[as.character(upgenes$row_num)]
    downgenesexp = expregion[as.character(downgenes$row_num)]
    uplimit = upperlimit[as.character(upgenes$row_num)]
    downlimit = lowerlimit[as.character(downgenes$row_num)]
    
    # Add disorder, region, mean, limit 
    upgenes = cbind(rep(disorder,nrow(upgenes)),rep(region,nrow(upgenes)),upgenes,rep("UP",nrow(upgenes)),upgenesexp,uplimit,rep(pvalue,nrow(upgenes)),rep(corrected,nrow(upgenes)))
    downgenes = cbind(rep(disorder,nrow(downgenes)),rep(region,nrow(downgenes)),downgenes,rep("DOWN",nrow(downgenes)),downgenesexp,downlimit,rep(pvalue,nrow(downgenes)),rep(corrected,nrow(downgenes)))
    colnames(upgenes) = c("disorder","region","row_num","gene_id","ensembl_gene_id","gene_symbol","entrez_id","upORdown","regionExpression","8PWDto37PCWmean","pvalue","corrected")
    colnames(downgenes) = c("disorder","region","row_num","gene_id","ensembl_gene_id","gene_symbol","entrez_id","upORdown","regionExpression","8PWDto37PCWmean","pvalue","corrected")
    tmp = rbind(upgenes,downgenes)

   # Finally, append the list to our result object
   results = rbind(results,tmp)
  }
}

# Save results to outfile!
save(results,file=resultfile)
