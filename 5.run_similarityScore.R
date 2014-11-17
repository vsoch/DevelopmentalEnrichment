# This is a submission script to run instances of creating a similarity matrix for clustering
# of gene -region sets (in 5.clusterRegionGenes.R)

# Directory to output Rda files (similarity vectors!)
outdir = "/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/reducedSimilarities"

setwd("/home/vsochat/SCRIPT/R/predictDisorder")

# Load the reduced feature matrix
load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/reduced_feature_matrix_51641.Rda")
data = temporal$matrix

jobsrun = 0
index = 1
while(jobsrun < nrow(data)){
  numjobs = as.numeric(system("squeue -u vsochat | wc -l",intern=TRUE)) - 2
  if (numjobs <= 4950){
    gene = rownames(data)[index]
    outfile = paste(outdir,"/",gene,"_simscores.Rda",sep="")
    # We need to print a job file for each term to run:
    if (!file.exists(outfile)) {
      jobby = paste("simscore_",gene,".job",sep="")
      sink(paste(".jobs/",jobby,sep=""))
      cat("#!/bin/bash\n")
      cat("#SBATCH --job-name=",jobby,"\n",sep="")  
      cat("#SBATCH --output=.out/",jobby,".out\n",sep="")  
      cat("#SBATCH --error=.out/",jobby,".err\n",sep="")  
      cat("#SBATCH --time=1-00:00\n",sep="")
      cat("#SBATCH --mem=12000\n",sep="")
      cat("Rscript /home/vsochat/SCRIPT/R/predictDisorder/5.similarityMatrix.R",gene,outfile,"\n")
      sink()
      # SUBMIT R SCRIPT TO RUN ON CLUSTER  
      system(paste("sbatch -p dpwall ",paste(".jobs/",jobby,sep="")))
      jobsrun = jobsrun + 1
      cat("Jobs run is",jobsrun,". Index is",index,"gene is",gene,"\n")
      index = index + 1
    } else {
      index = index + 1
    }
  }
}
