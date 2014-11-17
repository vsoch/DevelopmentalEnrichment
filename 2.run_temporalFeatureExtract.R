# This is a submission script to run instances of matching NeuroSynth Topic Maps to points
# in the Allen Brain Atlas

# Directory to save png images
imgdir = "/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/timeseries_images"

# Directory to output Rda files (features!)
outdir = "/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/timeseries_features"

# Load the gene timeseries
load("/scratch/PI/dpwall/DATA/ALLEN-BRAIN/DEVELOPING/genes_47808_timeseries.Rda")

for (f in 25001:length(allseries)){
  gene = names(allseries)[f]
  outfile = paste(outdir,"/",gene,"_temporalFeatures.Rda",sep="")
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
    cat("Rscript /home/vsochat/SCRIPT/R/predictDisorder/temporalFeatureExtract.R",gene,imgdir,outfile,"\n")
    sink()
    # SUBMIT R SCRIPT TO RUN ON CLUSTER  
    system(paste("sbatch -p dpwall ",paste(".jobs/",jobby,sep="")))
  }
}
