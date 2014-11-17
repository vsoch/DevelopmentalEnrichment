# This script will read in the reduced temporal features for ~50K gene probes, and parse into one matrix!  We hope that this matrix is much smaller than the original with ~700K rows

setwd("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/reducedSimilarities")
files = list.files(pattern="*_simscores.Rda")
load(files[1])

# Output folder
outdir = "/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/finalSimMatrices"

# We will parse rows into different matrices
scores = names(sim)
scores = scores[-which(scores %in% c("gene","rowidx"))]

# Create matrices to hold values
cat("Creating matrices...\n")
matrix.euc = array(dim=c(length(sim$binary),length(sim$binary)))
matrix.cos = array(dim=c(length(sim$binary),length(sim$binary)))
matrix.max = array(dim=c(length(sim$binary),length(sim$binary)))
matrix.man = array(dim=c(length(sim$binary),length(sim$binary)))
matrix.can = array(dim=c(length(sim$binary),length(sim$binary)))
matrix.bin = array(dim=c(length(sim$binary),length(sim$binary)))
matrix.min = array(dim=c(length(sim$binary),length(sim$binary)))
rownames(matrix.euc) = names(sim$binary)
colnames(matrix.euc) = names(sim$binary)
rownames(matrix.cos) = names(sim$binary)
colnames(matrix.cos) = names(sim$binary)
rownames(matrix.max) = names(sim$binary)
colnames(matrix.max) = names(sim$binary)
rownames(matrix.man) = names(sim$binary)
colnames(matrix.man) = names(sim$binary)
rownames(matrix.can) = names(sim$binary)
colnames(matrix.can) = names(sim$binary)
rownames(matrix.bin) = names(sim$binary)
colnames(matrix.bin) = names(sim$binary)
rownames(matrix.min) = names(sim$binary)
colnames(matrix.min) = names(sim$binary)

# Load each file and add to matrix
cat("Adding scores to matrix...\n")
for (f in 1:length(files)){
  load(files[f])
  cat("Adding scores in file",f,"\n")
  #matrix.cos[sim$gene,] = sim$cosine
  #matrix.cos[sim$gene,] = sim$cosine
  #matrix.max[sim$gene,] = sim$maximum
  #matrix.man[sim$gene,] = sim$manhattan
  matrix.can[sim$gene,] = sim$canberra
  #matrix.bin[sim$gene,] = sim$binary
  #matrix.min[sim$gene,] = sim$minkowski
}

cat("Saving matrices to file...\n")
#save(matrix.euc,file=paste(outdir,"/euclidean_similarity.Rda",sep=""))
#save(matrix.cos,file=paste(outdir,"/cosine_similarity.Rda",sep=""))
#save(matrix.max,file=paste(outdir,"/maximum_similarity.Rda",sep=""))
#save(matrix.man,file=paste(outdir,"/manhattan_similarity.Rda",sep=""))
save(matrix.can,file=paste(outdir,"/canberra_similarity.Rda",sep=""))

save(matrix.bin,file=paste(outdir,"/binary_similarity.Rda",sep=""))
save(matrix.min,file=paste(outdir,"/minkowski_similarity.Rda",sep=""))


