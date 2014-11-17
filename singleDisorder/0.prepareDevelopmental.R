# This script will parse the Allen Brain Atlas developmental data
# VSochat: September 27, 2014

# Let's first just get them into matrices

# EXOME --> protein coding portion of a gene ---------------
# EXON - GENES
setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/exome/exon/genes")
raw = read.csv("expression_matrix.csv",head=TRUE,sep=",")
# Data matrix is 17603 (genes, rows) by 493 (people,columns)

subjects = read.csv("columns_metadata.csv",head=TRUE,sep=",")
genes = read.csv("rows_metadata.csv",head=TRUE,sep=",")
exon.genes = list(columns=subjects,rows=genes,data=raw)
save(exon.genes,file="../../../exon-genes-17603x493.Rda")

# EXON - PROBES
setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/exome/exon/probes")
raw = read.csv("expression_matrix.csv",head=TRUE,sep=",")
# Data matrix is 17603 (genes, rows) by 493 (people,columns)

subjects = read.csv("columns_metadata.csv",head=TRUE,sep=",")
genes = read.csv("rows_metadata.csv",head=TRUE,sep=",")
exon.probes = list(columns=subjects,rows=genes,data=raw)
# Here we have 230694 probes, 493 brains
save(exon.probes,file="../../../exon-probes-230694x493.Rda")

# RNA GENES --> 
setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/exome/rna/genes")
raw = read.csv("expression_matrix.csv",head=TRUE,sep=",")
# Data matrix is 52375 (genes, rows) by 493 (people,columns)

subjects = read.csv("columns_metadata.csv",head=TRUE,sep=",")
genes = read.csv("rows_metadata.csv",head=TRUE,sep=",")
rna.genes = list(columns=subjects,rows=genes,data=raw)
save(rna.genes,file="../../../rna-genes-52375x525.Rda")

# RNA EXONS --> 
setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/exome/rna/exons")
raw = read.csv("expression_matrix.csv",head=TRUE,sep=",")
# Data matrix is 309222 (genes, rows) by 493 (people,columns)

subjects = read.csv("columns_metadata.csv",head=TRUE,sep=",")
genes = read.csv("rows_metadata.csv",head=TRUE,sep=",")
rna.exons = list(columns=subjects,rows=genes,data=raw)
save(rna.exons,file="../../../rna-exons-309222x525.Rda")


# MRNA ---------------------------------------------------------
setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/mrna")
# It looks like we have data for two subjects - 1109,1110

# lines 1 through 14 are meta data
# line 15 starts the data matrix
raw = read.csv("1109_methylation_beta_values.txt",head=TRUE,skip=14,sep="\t")
meta = readLines("1109_methylation_beta_values.txt")
meta = meta[1:14]
# 485577 by 87

tmp = c()
columns = c()
for (m in 1:length(meta)){
  line = meta[m]
  line = strsplit(gsub('"',"",line),"\t")[[1]]
  tmp = cbind(tmp,line[2:length(line)])
  columns = c(columns,line[1])
}

tmp = as.data.frame(tmp)
colnames(tmp) = columns
mrna.1109 = list(meta=tmp,data=raw)
save(mrna.1109,file="../../../mrna-1109-485577x87.Rda")

raw = read.csv("1110_methylation_beta_values.txt",head=TRUE,skip=14,sep="\t")
meta = readLines("1110_methylation_beta_values.txt")
meta = meta[1:14]
# 485577 by 92

tmp = c()
columns = c()
for (m in 1:length(meta)){
  line = meta[m]
  line = strsplit(gsub('"',"",line),"\t")[[1]]
  tmp = cbind(tmp,line[2:length(line)])
  columns = c(columns,line[1])
}

tmp = as.data.frame(tmp)
colnames(tmp) = columns
mrna.1110 = list(meta=tmp,data=raw)
save(mrna.1110,file="../../../mrna-1109-485577x92.Rda")

# Not sure what is in this file vs the others.
raw = read.csv("MicroRNA.csv",head=TRUE,sep=",")
structures = read.csv("structures.tab",head=TRUE,sep="\t")
specimens = read.csv("specimen.tab",head=TRUE,sep="\t")
mrna = list(data=raw,specimens=specimens,structures=structures)
save(mrna,file="../../../mrna-1861x217.Rda")

# TRACT ---------------------------------------------------------
# TO DO - need to figure out how to turn this into features
setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/tract")


# LMD ---------------------------------------------------------
# TO DO - need to figure out how to turn this into features
setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/LMD")
files = list.files()

# let's get everything in Rda'
for (f in 1:length(files)){
  cat("Processing",f,"of",length(files),"\n")
  file = files[f]
  subid = strsplit(file,"_")[[1]][1]
  structure = strsplit(file,"[.]")[[1]][3]
  structure = paste(strsplit(structure,"_")[[1]][-1],collapse="_")
  structure = gsub("_LMD","",structure)
  # Read in the file - data starts at 10th line, 1-10 are meta
  raw = read.csv(file,head=TRUE,sep="\t",skip=9)
  meta = readLines(file)
  meta = meta[1:9]
  # Parse the meta portion
  feparams = strsplit(meta[2],"\t")[[1]][-1]
  fe = strsplit(meta[3],"\t")[[1]][-1]
  statparams = strsplit(meta[6],"\t")[[1]][-1]
  stats = strsplit(meta[7],"\t")[[1]][-1]
  names(stats) = statparams
  names(fe) = feparams
  # Save to file
  person = paste(subid,structure,sep="[.]")
  statsall[[person]] = stats
  featuresall[[person]] = fe
  lmd = list(data=raw,subid=subid,file=file,structure=structure,stats=stats,features=fe)
  save(lmd,file=paste(person,".Rda",sep=""))
}
