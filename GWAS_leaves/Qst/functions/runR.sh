#!/bin/bash

##SBATCH -t unlimited (for tweedie)
##SBATCH -p nolimit (for tweedie)
##SBATCH --mem=10GB (for tweedie)
#SBATCH -t 48:00:00
#SBATCH -p core
#SBATCH -n 1
#SBATCH --mail-type FAIL
#SBATCH --mail-user=sara.westman@umu.se
#SBATCH --account u2015011

module load R/4.1.1
      
Rscript $@