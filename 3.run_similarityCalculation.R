# This is a submission script to run instances of matching NeuroSynth Topic Maps to points
# in the Allen Brain Atlas

setwd("/home/vsochat/SCRIPT/R/predictDisorder")

# Directory to output Rda files (similarity vectors!)
outdir = "/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/geneRegionSimilarity"

# Load the gene feature matrix
load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/timeseries_featurematrix.Rda")

# Load the gene list
load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/timeseries_genes_list.Rda")

jobsrun = 0
index = 1
while(jobsrun < length(genes)){
  numjobs = as.numeric(system("squeue -u vsochat | wc -l",intern=TRUE)) - 2
  if (numjobs <= 4950){
    gene = genes[index]
    outfile = paste(outdir,"/",gene,"_reduced.Rda",sep="")
    # We need to print a job file for each term to run:
    if (!file.exists(outfile)) {
      jobby = paste(gene,".job",sep="")
      sink(paste(".jobs/",jobby,sep=""))
      cat("#!/bin/bash\n")
      cat("#SBATCH --job-name=",jobby,"\n",sep="")  
      cat("#SBATCH --output=.out/",jobby,".out\n",sep="")  
      cat("#SBATCH --error=.out/",jobby,".err\n",sep="")  
      cat("#SBATCH --time=1-00:00\n",sep="")
      cat("#SBATCH --mem=12000\n",sep="")
      cat("Rscript /home/vsochat/SCRIPT/R/predictDisorder/3.similarityCalculation.R",gene,outfile,"\n")
      sink()
      # SUBMIT R SCRIPT TO RUN ON CLUSTER  
      system(paste("sbatch -p dpwall ",paste(".jobs/",jobby,sep="")))
      jobsrun = jobsrun + 1
      cat("Jobs run is",jobsrun,". Index is",index,"\n")
      index = index + 1
    } 
  }
}