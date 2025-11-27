#######################################################################
# runfastJT 
#######################################################################
# Load packages 
suppressPackageStartupMessages({
  library(dplyr)
  library(fastJT)
  library(vroom)
  library(qqman)
  library(tidyverse)
  library(tibble)
  library(utils)
  library(stringr)
  })

# Set parameters 
out_dir <- "/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/fastJT/res_tweedie"
file_name <- "SPG_leaves_tweedie"
pheno <- "/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/Metabolite_data/BLUPs/SPGs_tweedie/SPGs_tweedie_leaves_BLUP.txt"

snpmat <-"/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_ML/data/SNPs/SNPmat_for_fastJT/SwAsp_SNPmat_fastJT.rds"
swasp99_genotypes <- read.table("/mnt/picea/projects/aspseq/nstreet/swasp/gwas/data/PLINK/SwAsp_AfterBatchRemoval_Het.MAF.HWE.fam", quote="\"", comment.char="") %>% dplyr::select(V2)

# Print info
message("Arguments provided")
message(paste("\tOutput directory:", out_dir))
message(paste("\tPhenotype data:", pheno))
message(paste("\tFile prefix:", file_name))
message(paste("\tSNP matrix:", snpmat))

# Load data 
geno_feature <- readRDS(snpmat)
pheno_dat <- read.delim(pheno)
  
# Prep metabolite data by selecting only those genotypes that are in the SNP data 
colnames(swasp99_genotypes) <- "Genotype"

pheno_dat$Genotype <- as.integer(pheno_dat$Genotype)
pheno_dat <- left_join(swasp99_genotypes, pheno_dat, by = "Genotype")
pheno_dat <- column_to_rownames(pheno_dat, var = "Genotype")
pheno_dat <- pheno_dat %>% as.matrix()

# Run fastJT 
fastJT_res <- fastJT(Y = pheno_dat, X = geno_feature, outTopN = NA)

# Run fastJT select
fastJT_res_sel <- fastJT.select(Y = pheno_dat, X = geno_feature, cvMesh = NULL, kFold = 10,
                                selCrit = NULL, outTopN = NA, numThreads = 10)

# Save data 
save(fastJT_res, fastJT_res_sel, file = file.path(out_dir, paste0(file_name, ".RData")))

#######################################################################
# Save session info
#######################################################################

writeLines(capture.output(sessionInfo()), con = file.path(out_dir, paste("sessionInfo",Sys.Date() ,".txt", sep="")))
message("Done")
