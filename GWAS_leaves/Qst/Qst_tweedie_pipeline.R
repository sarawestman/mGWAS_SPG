suppressPackageStartupMessages({
  library(dplyr)
  library(cplm)
})
source("/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/scripts/R/Pop_genetics/Qst/tweedie/generate_Qst_tweedie.R")

#######################################################################
# Bash loop input arguments 
#######################################################################
args <- commandArgs(trailingOnly = TRUE)
cat(args, sep = "\n")

# test if there is two argument: if not, return an error
if (length(args) != 3) {
  stop("Three argument must be supplied", call.=FALSE)
} 

# Direct the arguments
tweedie_colnr <- args[1] %>% as.integer()
trait_name <- args[2]
out_dir <- args[3]

#######################################################################
# Set arguments 
#######################################################################
# Data 
input_data <- "/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/Pop_genetics/Qst/tweedie_extras/SPGs_buds_Outlier_removed_tweedie_extras_data.3.reps_pop_lat.txt"

# generate_Qst_tweedie arguments - specify column number. See generate_Qst_tweedie.R function for more info about each parameter.   
random_effect = c(2,4) 
fixed_effect = 3
start_column = 5
col_convert = TRUE
fixef_contrasts = TRUE
bootstr_R = 100 # 100 for Coumarin_glycoside_1__Fraxin_X_3.2144_369.0824 and Feruloylquinic_acid_hexoside_X_5.3761_529.1547, 1000 for the rest

# Due to loop - Adjust column numbers to match the input_data column order (before the first trait column)
raw_data <- readr::read_tsv(input_data)[,c(1:4,tweedie_colnr)]

#######################################################################
# Generate Qst
#######################################################################

# Print arguments 
message("Arguments provided")
message(paste("\tInput data:", input_data))
message(paste("\tOutput directory:", out_dir))
message(paste("\tRandom effect:", colnames(raw_data)[random_effect]))
message(paste("\tFixed effect:", colnames(raw_data)[fixed_effect]))
message(paste("\tTrait column:", colnames(raw_data)[start_column]))
message(paste("\tBootstrap posterior sample size:", bootstr_R))

# Create termlabel formula
termlabs <- c(paste0("(1|", colnames(raw_data)[random_effect], ")"))
if (length(random_effect) > 0 & !is.null(fixed_effect)) { # with fixed effect
  termlabs <- c(colnames(raw_data)[fixed_effect], termlabs)
}

# Generate Qst
results <-
  generate_Qst_tweedie(dat = raw_data, random_effect = random_effect, start_column = start_column, termlabs = termlabs,
                       fixed_effect = fixed_effect, col_convert = col_convert, 
                       fixef_contrasts = fixef_contrasts, bootstr_R = bootstr_R)

# Save data
readr::write_tsv(results, file = file.path(out_dir, paste0(trait_name, "_Qst_table.txt"))) 

#######################################################################
# Save session info
#######################################################################
writeLines(capture.output(sessionInfo()), con = file.path(out_dir, paste("sessionInfo",Sys.Date() ,".txt", sep="")))
message("Done")