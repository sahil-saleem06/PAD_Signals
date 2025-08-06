#!/bin/bash

#BSUB -J GWAS_Inspector_pad_all_gr38
#BSUB -n 1
#BSUB -R rusage[mem=64000]
#BSUB -q voltron_normal
#BSUB -W 48:00
#BSUB -o GWAS_Inspector_pad_all_gr38.out
#BSUB -e GWAS_Inspector_pad_all_gr38.err

# set R library directory
export R_LIBS_USER=$HOME/R/rocker-rstudio/bioconductor-tidyverse_3.17

# ensure singularity can see the correct folders
export SINGULARITY_BIND="/lsf:/lsf, /project/:/project/, /appl/:/appl/, /lsf/:/lsf/, /scratch/:/scratch, /static:/static"

set -e

module load singularity

singularity exec /project/voltron/rstudio/containers/bioconductor-tidyverse_3.17.sif Rscript /project/damrauer_scratch/Users/saleemsa/PAD_Signals/GWASInspectorFiles/GWAS.ins.R