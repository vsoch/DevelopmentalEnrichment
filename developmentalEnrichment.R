## DEVELOPMENTAL BRAIN
# Is there a timepoint in development that has enrichment for disorder genes in the up and down?

args <- commandArgs(TRUE)
disorder = args[1]
resultfile = args[2]

load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/rna-genes-52375x525.Rda")
# Here are disorder genes
load("/scratch/users/vsochat/DATA/ALLEN/UberDisorder/Uber5_NonZeroGenes.Rda")

# Let's look at the mean and standard deviation for all, in utero, child, teenager, adult
#utero = c("8 pcw","9 pcw","12 pcw","13 pcw","16 pcw","17 pcw","19 pcw","21 pcw","24 pcw","25 pcw","26 pcw","35 pcw","37 pcw")
#child = c("4 mos","10 mos","1 yrs","2 yrs","3 yrs","4 yrs","8 yrs","11 yrs","13 yrs","15 yrs") 
#adult = c("18 yrs","19 yrs","21 yrs","23 yrs","30 yrs","36 yrs","37 yrs","40 yrs")

#uidx = rna.genes$columns[which(rna.genes$columns$age %in% utero),]
#subset = rna.genes$data[,uidx$column_num]
#mean(as.numeric(subset[1,]))
#sd(as.numeric(subset[1,]))

#uidx = rna.genes$columns[which(rna.genes$columns$age %in% child),]
#subset = rna.genes$data[,uidx$column_num]
#mean(as.numeric(subset[1,]))
#sd(as.numeric(subset[1,]))

#uidx = rna.genes$columns[which(rna.genes$columns$age %in% adult),]
#subset = rna.genes$data[,uidx$column_num]
#mean(as.numeric(subset[1,]))
#sd(as.numeric(subset[1,]))

# We are going to define means/sd based on development - the distribution is higher
# and more varied than later in life, and we are interested if genes are more expressed
# at a certain timepoint

# This is the time period we will filter to
utero = c("8 pcw","9 pcw","12 pcw","13 pcw","16 pcw","17 pcw","19 pcw","21 pcw","24 pcw","25 pcw","26 pcw","35 pcw","37 pcw")
uidx = rna.genes$columns[which(rna.genes$columns$age %in% utero),]
people = rna.genes$columns[which(rna.genes$columns$age %in% utero),]
probes = rna.genes$rows
probes$gene_symbol =  as.character(probes$gene_symbol)

data = rna.genes$data[,uidx$column_num]

# Let's calculate the mean for each gene (rows)
genemeans = rowMeans(data)
genesd = apply(data,1,sd)
upperlimit = genemeans + 3*genesd
lowerlimit = genemeans - 3*genesd

setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING")
timepoints = utero

# ANALYSIS 2  
# This first analysis will look at the gene across all regions, for each timepoint
# The mean being compared to is the mean and sd across ALL DEVELOPMENT

fisher = c()
for (g in 1:length(genes)){
  disorder = names(genes[g])
  cat("Processing disorder",disorder,"\n")
  gen = names(genes[[disorder]])

  disordergenes = probes$gene_symbol[which(probes$gene_symbol %in% gen)]
  nondisordergenes = probes$gene_symbol[-which(probes$gene_symbol %in% gen)]
  
  for (t in 1:length(timepoints)){
    cat("Processing timepoint",t,"\n")
    # Get data for that timepoint
    exp = data[,which(people$age %in% timepoints[t])]
    ppl = people[which(people$age %in% timepoints[t]),]

    # Find genes above/below 3SD of mean
    if (!is.vector(exp)){
      weird = probes[c(which(rowMeans(exp) >= upperlimit),which(rowMeans(exp) >= upperlimit)),]
      notweird = probes[intersect(which(rowMeans(exp) <= upperlimit),which(rowMeans(exp) >= lowerlimit)),]
    } else{
      weird = probes[c(which(exp >= upperlimit),which(exp >= upperlimit)),]
      notweird = probes[intersect(which(exp <= upperlimit),which(exp >= lowerlimit)),]
    }  
    # Make a 2x2 table for the regions/disorder
    # High or Low Expression and disorder gene
    a = length(which(weird$gene_symbol %in% disordergenes))
    # Not high or low expression and disorder gene
    b = length(which(notweird$gene_symbol %in% disordergenes))
    # High or Low Expression and non-disorder gene
    c = length(which(weird$gene_symbol %in% nondisordergenes))
    # Not high or low expression and non-disorder gene
    d = length(which(notweird$gene_symbol %in% nondisordergenes))

    tabley = matrix(c(a,c,b,d),nrow=2,dimnames=list(c("DisorderGenes","~DisorderGenes"),c("WeirdExpression","~WeirdExpression")))

    ft = fisher.test(tabley, conf.level = 0.95)
    fisher = rbind(fisher,c(disorder,timepoints[t],ft$p.value))
  }
}

fisher = as.data.frame(fisher,stringsAsFactors=FALSE)
colnames(fisher) = c("DISORDER","SAMPLEID","PVALUE")
fisher$PVALUE = as.numeric(fisher$PVALUE)

# Let's FDR correct'
FDR = p.adjust(fisher$PVALUE)
# Save fisher result to outfile
fisher = cbind(fisher,FDR)
save(fisher,file="developingBrainTimepointsFisher.Rda")

# Find significant results!
which(fisher$FDR<=0.05)

# ANALYSIS 3
# This next analysis will look on a regional basis - is the expression in a region 
# of the brain different than other parts of the brain for a particular timepoint?
# This is the same as the first analysis, but we have multiple timepoints
# We would want to see some kind of clustering of regions at a particular timepoint

# We are again going to compare to the means across ALL DEVELOPMENT
genemeans = rowMeans(data)
genesd = apply(data,1,sd)
upperlimit = genemeans + 3*genesd
lowerlimit = genemeans - 3*genesd

fisher = c()
# For each disorder
for (g in 1:length(genes)){
  disorder = names(genes[g])
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

      # Find genes above/below 3SD of mean
      #if (!is.vector(exp)){
      #  weird = probes[c(which(rowMeans(exp) >= upperlimit),which(rowMeans(exp) >= upperlimit)),]
      #  notweird = probes[intersect(which(rowMeans(exp) <= upperlimit),which(rowMeans(exp) >= lowerlimit)),]
      #} else{
      weird = probes[c(which(expregion >= upperlimit),which(expregion >= upperlimit)),]
      notweird = probes[intersect(which(expregion <= upperlimit),which(expregion >= lowerlimit)),]
      #}  
      # Make a 2x2 table for the regions/disorder
      # High or Low Expression and disorder gene
      a = length(which(weird$gene_symbol %in% disordergenes))
      # Not high or low expression and disorder gene
      b = length(which(notweird$gene_symbol %in% disordergenes))
      # High or Low Expression and non-disorder gene
      c = length(which(weird$gene_symbol %in% nondisordergenes))
      # Not high or low expression and non-disorder gene
      d = length(which(notweird$gene_symbol %in% nondisordergenes))

      tabley = matrix(c(a,c,b,d),nrow=2,dimnames=list(c("DisorderGenes","~DisorderGenes"),c("WeirdExpression","~WeirdExpression")))

      ft = fisher.test(tabley, conf.level = 0.95)
      fisher = rbind(fisher,c(disorder,timepoints[t],region,ft$p.value))
    }
  }
}
fisher = as.data.frame(fisher,stringsAsFactors=FALSE)
colnames(fisher) = c("DISORDER","SAMPLEID","PVALUE")
fisher$PVALUE = as.numeric(fisher$PVALUE)
# Save fisher result to outfile
save(fisher,file=resultfile)


}


# Can I instead look at change between stages of development?

