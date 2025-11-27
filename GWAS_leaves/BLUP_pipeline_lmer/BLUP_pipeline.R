suppressPackageStartupMessages({
  library(dplyr)
})
source("/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/scripts/R/Pre_GWAS/BLUP_pipeline/generate_BLUP.R")
source("/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/scripts/R/Pre_GWAS/BLUP_pipeline/plot_ranef_title.R")

#######################################################################
# Set arguments 
#######################################################################
input_data <- "/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/Metabolite_data/raw_data/s12.SPGs.final.Dec2019_noNA_meta_data.txt"
out_dir <- "/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/Metabolite_data/BLUPs/SPGs"
data_name <- "SPGs_leaves" # Filename prefix
line_name <- "Clone" # Has to match line column name in input_data. E.g., "Genotype", "Clone". Required for boxplots. 

# generate_BLUP arguments - specify column number. See generate_BLUP.R function for more info about each parameter.   
sample_column = 1
random_effect = 7
fixed_effect = 9
start_column = 10
col_convert = TRUE
fixef_contrasts = TRUE

#######################################################################
# Generate BLUPs
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

# Generate BLUP
results <-
  generate_BLUP(dat = raw_data, sample_column = sample_column, random_effect = random_effect, start_column = start_column, 
                fixed_effect = fixed_effect, col_convert = col_convert, fixef_contrasts = fixef_contrasts)

if (is.list(results)) {
  readr::write_tsv(results$BLUP, file = file.path(out_dir, paste0(data_name, "_BLUP.txt")))
  saveRDS(results$Outliers_residuals, file = file.path(stats_out_dir, paste0(data_name, "_Outliers_residuals.rds")))
  readr::write_tsv(results$Outlier_data, file = file.path(stats_out_dir, paste0(data_name, "_Outlier_data.txt")))
  readr::write_tsv(results$Outlier_removed_data, file = file.path(stats_out_dir, paste0(data_name, "_Outlier_removed_data.txt")))
  readr::write_tsv(results$Shapiro_raw_data, file = file.path(stats_out_dir, paste0(data_name, "_raw_shapiro.txt")))
  readr::write_tsv(results$Shapiro_BLUP, file = file.path(stats_out_dir, paste0(data_name, "_BLUP_shapiro.txt")))
  readr::write_tsv(results$Norm_raw_dat, file = file.path(stats_out_dir, paste0(data_name, "_Outlier_rm_norm_trans_data.txt")))
  readr::write_tsv(results$Transformed_traits, file = file.path(stats_out_dir, paste0(data_name, "_transformed_traits.txt")))
  readr::write_tsv(results$Non_transformed_traits, file = file.path(stats_out_dir, paste0(data_name, "_non_transformed_traits.txt")))
  if(exists('Removed_traits', where = results)) {
    readr::write_tsv(results$Removed_traits, file = file.path(stats_out_dir, paste0(data_name, "_removed_traits.txt")))
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
  ggplot2::ggsave(file.path(stats_plot_dir, "BLUP_histogram.pdf"), ml, width = 20, height = 15)
}

# Boxplots - color outliers
if (is.list(results)) {
  outlier_list <- results$Outliers_residuals
  trait_names <- outlier_list[lapply(outlier_list,length)>0] %>% names()
  raw_data$my_color <- "NULL"
  # Convert lines to factor type
  raw_data[[line_name]] <- as.factor(raw_data[[line_name]])
  # Generate plots 
  cbPalette = c("darkcyan", "coral1")
  plot_lst <- lapply(trait_names, function(t) {
    my_rows <- outlier_list[[t]]
    for(i in 1:nrow(raw_data)){
      if(i %in% my_rows){
        raw_data[i, "my_color"] <- "rm_outliers"
      } else {
        raw_data[i, "my_color"] <- "kept_samples"
      }
    }
    p <- ggplot2::ggplot(raw_data, ggplot2::aes(x = .data[[line_name]], y = .data[[t]])) +
      ggplot2::geom_boxplot() + 
      ggplot2::geom_jitter(ggplot2::aes(color = `my_color`)) + 
      ggplot2::scale_fill_manual(values = cbPalette) +
      ggplot2::scale_color_manual(values = cbPalette) +
      ggplot2::labs(title = t) 
    return(p)
  })  
  message("Creating boxplots")
  ml <- gridExtra::marrangeGrob(plot_lst, nrow = 1, ncol = 1)
  ggplot2::ggsave(file.path(stats_plot_dir, "Boxplots_residual_outlier_rm.pdf"), ml, width = 25, height = 12)
}

# Model parameters - nontransformed traits 
if (is.list(results)) {
  lmms_all <- results$lme_data
  # Extract nontransformed trait names 
  non_trans_trait <- names(lmms_all)[names(lmms_all) %in% results$Non_transformed_traits$non_trans_trait]
  # Generate plots 
  plot_lst <- lapply(non_trans_trait, function(t) {
    # creates a residual quantile plot for the error term
    my_plot <- redres::plot_resqq(lmms_all[[t]]) + ggplot2::ggtitle(t) 
  })
  plot_lst2 <- lapply(non_trans_trait, function(t) {
    # creates normal quantile plots for each random effect
    #my_plot <- redres::plot_ranef(lmms_all[[t]]) 
    plot_ranef_title(lmms_all[[t]], trait_title = t)
  })
  plot_lst3 <- lapply(non_trans_trait, function(t) {
    # creates a plot of the conditional studentized residuals versus the fitted values
    my_plot <- redres::plot_redres(lmms_all[[t]], type = "std_cond") + ggplot2::ggtitle(t)
  })
  message("Creating model assumption plots (non-transformed)")
  ml1 <- gridExtra::marrangeGrob(plot_lst, nrow = 5, ncol = 5)
  ml2 <- gridExtra::marrangeGrob(plot_lst2, nrow = 5, ncol = 5)
  ml3 <- gridExtra::marrangeGrob(plot_lst3, nrow = 5, ncol = 5)
  ggplot2::ggsave(file.path(stats_plot_dir, "Residual_qqplots_nontransformed_traits.pdf"), ml1, width = 25, height = 12)
  ggplot2::ggsave(file.path(stats_plot_dir, "Random_effect_qqplots_nontransformed_traits.pdf"), ml2, width = 25, height = 12)
  ggplot2::ggsave(file.path(stats_plot_dir, "Homoscedasticity_nontransformed_traits.pdf"), ml3, width = 25, height = 12)
  }

# Model parameters - transformed traits 
if (is.list(results)) {
  # Extract transformed trait names 
  trans_trait <- names(lmms_all)[names(lmms_all) %in% results$Transformed_traits$trans_trait]
  # Generate plots 
  plot_lst <- lapply(trans_trait, function(t) {
    # creates a residual quantile plot for the error term
    my_plot <- redres::plot_resqq(lmms_all[[t]]) + ggplot2::ggtitle(t) 
  })
  plot_lst2 <- lapply(trans_trait, function(t) {
    # creates normal quantile plots for each random effect
    my_plot <- plot_ranef_title(lmms_all[[t]], trait_title = t)
  })
  plot_lst3 <- lapply(trans_trait, function(t) {
    # creates a plot of the conditional studentized residuals versus the fitted values
    my_plot <- redres::plot_redres(lmms_all[[t]], type = "std_cond") + ggplot2::ggtitle(t)
  })
  message("Creating model assumption plots (transformed traits)")
  ml1 <- gridExtra::marrangeGrob(plot_lst, nrow = 5, ncol = 5)
  ml2 <- gridExtra::marrangeGrob(plot_lst2, nrow = 5, ncol = 5)
  ml3 <- gridExtra::marrangeGrob(plot_lst3, nrow = 5, ncol = 5)
  ggplot2::ggsave(file.path(stats_plot_dir, "Residual_qqplots_transformed_traits.pdf"), ml1, width = 25, height = 12)
  ggplot2::ggsave(file.path(stats_plot_dir, "Random_effect_qqplots_transformed_traits.pdf"), ml2, width = 25, height = 12)
  ggplot2::ggsave(file.path(stats_plot_dir, "Homoscedasticity_transformed_traits.pdf"), ml3, width = 25, height = 12)
}

#######################################################################
# Save session info
#######################################################################

writeLines(capture.output(sessionInfo()), con = file.path(stats_out_dir, paste("sessionInfo",Sys.Date() ,".txt", sep="")))
message("Done")