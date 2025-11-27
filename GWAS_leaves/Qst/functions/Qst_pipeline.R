suppressPackageStartupMessages({
  library(dplyr)
})
source("/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/scripts/R/Pop_genetics/Qst/lmer/function/generate_Qst.R")
source("/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/scripts/R/Pop_genetics/Qst/lmer/function/plot_ranef_title.R")

#######################################################################
# Set arguments 
#######################################################################
input_data <- "/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/Pop_genetics/Qst/lmer/SPGs_buds_Outlier_removed_lmer_data.3.reps_pop_lat.txt" # OBS! To match the GWAS input I have normalized the same traits that were normalized for the GWAS  
out_dir <- "/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/Pop_genetics/Qst/lmer"
data_name <- "SPGs_leaves_lmer" # Filename prefix

# generate_Qst arguments - specify column number. See generate_Qst.R function for more info about each parameter.   
random_effect = c(2,4) 
fixed_effect = 3
start_column = 5
col_convert = TRUE
fixef_contrasts = TRUE
bootstr_R = 3 #just for testing, should be 1000

#######################################################################
# Generate Qst
#######################################################################
# Create stats folder and subfolder
stats_out_dir <- file.path(out_dir, "stats")
dir.create(stats_out_dir, recursive = TRUE, showWarnings = FALSE)
stats_plot_dir <- file.path(stats_out_dir, "plots")
dir.create(stats_plot_dir, recursive = TRUE, showWarnings = FALSE)

# Read file 
raw_data <- readr::read_tsv(input_data)

# Print arguments 
message("Arguments provided")
message(paste("\tInput data:", input_data))
message(paste("\tOutput directory:", out_dir))
message(paste("\tRandom effect:", colnames(raw_data)[random_effect]))
message(paste("\tFixed effect:", colnames(raw_data)[fixed_effect]))
message(paste("\tTrait start column:", colnames(raw_data)[start_column]))
message(paste("\tBootstrap posterior sample size:", bootstr_R))

# Generate Qst
results <-
  generate_Qst(dat = raw_data, random_effect = random_effect, start_column = start_column, 
               fixed_effect = fixed_effect, col_convert = col_convert, fixef_contrasts = fixef_contrasts, bootstr_R = bootstr_R)

if (is.list(results)) {
  readr::write_tsv(results$BLUP, file = file.path(stats_out_dir, paste0(data_name, "_Qst_BLUP.txt")))
  readr::write_tsv(results$Shapiro_BLUP, file = file.path(stats_out_dir, paste0(data_name, "Qst_BLUP_shapiro.txt")))
  readr::write_tsv(results$Qst_table, file = file.path(out_dir, paste0(data_name, "_Qst_table.txt"))) 
  if(exists('Removed_traits', where = results)) {
    readr::write_tsv(results$Removed_traits, file = file.path(stats_out_dir, paste0(data_name, "_Qst_removed_traits.txt")))
  }
}

#######################################################################
# Generate plots
#######################################################################
# Histogram BLUPs
if (is.list(results)) {
  BLUP <- results$BLUP
  # Generate plots 
  plot_lst <- lapply(names(BLUP[-1]), function(t) {
    hist_plot <- ggplot2::ggplot(BLUP, ggplot2::aes(.data[[t]])) + 
      ggplot2::geom_histogram(bins = 20) +
      ggplot2::ggtitle(t)
  })
  message("Creating BLUP histogram plots")
  ml <- gridExtra::marrangeGrob(plot_lst, nrow = 6, ncol = 6)
  ggplot2::ggsave(file.path(stats_plot_dir, "Qst_BLUP_histogram.pdf"), ml, width = 20, height = 15)
}

# Model parameters
if (is.list(results)) {
  lmms_all <- results$lme_data
  plot_lst <- lapply(names(lmms_all), function(t) {
    # creates a residual quantile plot for the error term
    my_plot <- redres::plot_resqq(lmms_all[[t]]) + ggplot2::ggtitle(t) 
  })
  plot_lst2 <- lapply(names(lmms_all), function(t) {
    # creates normal quantile plots for each random effect
    plot_ranef_title(lmms_all[[t]], trait_title = t)
  })
  plot_lst3 <- lapply(names(lmms_all), function(t) {
    # creates a plot of the conditional studentized residuals versus the fitted values
    my_plot <- redres::plot_redres(lmms_all[[t]], type = "std_cond") + ggplot2::ggtitle(t)
  })
  message("Creating model assumption plots")
  ml1 <- gridExtra::marrangeGrob(plot_lst, nrow = 5, ncol = 5)
  ml2 <- gridExtra::marrangeGrob(plot_lst2, nrow = 2, ncol = 2)
  ml3 <- gridExtra::marrangeGrob(plot_lst3, nrow = 5, ncol = 5)
  ggplot2::ggsave(file.path(stats_plot_dir,"Qst_Residual_qqplots.pdf"), ml1, width = 25, height = 12)
  ggplot2::ggsave(file.path(stats_plot_dir,"Qst_Random_effect_qqplots.pdf"), ml2, width = 25, height = 12)
  ggplot2::ggsave(file.path(stats_plot_dir,"Qst_Homoscedasticity.pdf"), ml3, width = 25, height = 12)
}

#######################################################################
# Save session info
#######################################################################
writeLines(capture.output(sessionInfo()), con = file.path(stats_out_dir, paste("sessionInfo",Sys.Date() ,".txt", sep="")))
message("Done")