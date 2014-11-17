# Set up script

setwd("/scratch/users/vsochat/DATA/ALLEN/UberDisorder")
load("Uber5_GenesGr1.Rda")
outdir = ("/scratch/users/vsochat/DATA/ALLEN/UberDisorder/developmentRegions1/")
setwd("/home/vsochat/SCRIPT/R/predictDisorder")

# Submission down here
for (f in 1:length(genes)){
  disorder = names(genes[f])
  # Output file with fisher test results
  outfile = paste(outdir,disorder,"_developmentRegionEnrichment.Rda",sep="")
  if (!file.exists(outfile)) {
    cat(disorder,"\n")
    jobby = paste(disorder,"dev.job",sep="")
    sink(paste(".jobs/",jobby,sep=""))
    cat("#!/bin/bash\n")
    cat("#SBATCH --job-name=",jobby,"\n",sep="")  
    cat("#SBATCH --output=.out/",jobby,".out\n",sep="")  
    cat("#SBATCH --error=.out/",jobby,".err\n",sep="")  
    cat("#SBATCH --time=0-03:00\n",sep="")
    cat("#SBATCH --mem=8000\n",sep="")
    cat("Rscript /home/vsochat/SCRIPT/R/predictDisorder/developmentEnrichment.R",disorder,outfile,"\n")
    sink()
  
    # SUBMIT R SCRIPT TO RUN ON CLUSTER  
    system(paste("sbatch",paste(".jobs/",jobby,sep="")))
  }
}


