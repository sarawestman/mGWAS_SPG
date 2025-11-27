suppressPackageStartupMessages({
  library(dplyr)
  library(pryr)
})
source("/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/scripts/R/Pre_GWAS/BLUP_tweedie_pipeline/generate_tweedie_BLUP.R")

#######################################################################
# Set arguments 
#######################################################################
input_data <- "/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/Metabolite_data/raw_data/s12.SPGs.final.Dec2019_noNA_meta_data.txt"
out_dir <- "/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/Metabolite_data/BLUPs/SPGs_tweedie"
data_name <- "SPGs_tweedie_leaves" # Filename prefix
line_name <- "Clone" # Has to match line column name in input_data. E.g., "Genotype", "Clone". Required for boxplots. 
trait_colnr <- read.delim("/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/Metabolite_data/BLUPs/SPGs_tweedie/Model_evaluation/Tweedie_trait_column_nr.tsv") %>% dplyr::pull(tweedie_colnr)

# generate_tweedie_BLUP arguments - specify column number. See generate_tweedie_BLUP.R function for more info about each parameter.   
random_effect = 7
fixed_effect = 9
trait_columns = trait_colnr
col_convert = TRUE
fixef_contrasts = TRUE
nsim = 10000

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
message(paste("\tNr of residual simulations:", nsim))
message(paste("\tTrait:", colnames(raw_data)[trait_columns], "\n"))

# Generate BLUP
results <-
  generate_tweedie_BLUP(dat = raw_data, random_effect = random_effect, trait_columns = trait_columns, nsim = nsim,
                fixed_effect = fixed_effect, col_convert = col_convert, fixef_contrasts = fixef_contrasts)

if (is.list(results)) {
  readr::write_tsv(results$BLUP, file = file.path(out_dir, paste0(data_name, "_BLUP.txt")))
  saveRDS(results$Outliers_residuals, file = file.path(stats_out_dir, paste0(data_name, "_Outliers_residuals.rds")))
  readr::write_tsv(results$Outlier_data, file = file.path(stats_out_dir, paste0(data_name, "_Outlier_data.txt")))
  readr::write_tsv(results$Outlier_removed_data, file = file.path(stats_out_dir, paste0(data_name, "_Outlier_removed_data.txt")))
  readr::write_tsv(results$testResidual_stats, file = file.path(stats_out_dir, paste0(data_name, "_testResidual_stats.txt")))
  readr::write_tsv(results$Shapiro_BLUP, file = file.path(stats_out_dir, paste0(data_name, "_BLUP_shapiro.txt")))
  readr::write_tsv(results$Zero_inflated_traits, file = file.path(stats_out_dir, paste0(data_name, "_zero_inflated_traits.txt")))
  readr::write_tsv(results$Not_zero_inflated_traits, file = file.path(stats_out_dir, paste0(data_name, "_not_zero_inflated_traits.txt")))
  readr::write_tsv(results$NotNormDist_BLUP_traits, file = file.path(stats_out_dir, paste0(data_name, "_notNormDist_BLUP_traits.txt")))
  readr::write_tsv(results$NormDist_BLUP_traits, file = file.path(stats_out_dir, paste0(data_name, "_normDist_BLUP_traits.txt")))
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
  ggplot2::ggsave(file.path(stats_plot_dir, "Tweedie_BLUP_histogram.pdf"), ml, width = 20, height = 15)
}

# Boxplots - color outliers
if (is.list(results)) {
  outlier_list <- results$Outliers_residuals
  trait_names <- outlier_list[lapply(outlier_list, length)>0] %>% names()
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

# Trait median/mean vs BLUPs 
if (is.list(results)) {
  BLUP <- results$BLUP
  dat <- results$Outlier_removed_data

  # Generate plots 
  plot_lst <- lapply(names(BLUP[-1]), function(t) {
    dat$trait <- dat[[t]]
    BLUP$blup_trait <- BLUP[[t]]
    
    # Prep median data 
    dat_median <- dat %>% dplyr::select(!!line_name, trait) %>% dplyr::group_by(!!!syms(line_name)) %>% dplyr::summarise(median(trait)) %>% data.table::setnames("median(trait)", "Trait_median")
    names(dat_median) <- c("Genotype", "Trait_median")
    dat_median$Genotype <- as.character(dat_median$Genotype)
    # Prep mean data 
    dat_mean <- dat %>% dplyr::select(!!line_name, trait) %>% dplyr::group_by(!!!syms(line_name)) %>% dplyr::summarise(mean(trait)) %>% data.table::setnames("mean(trait)", "Trait_mean")
    names(dat_mean) <- c("Genotype", "Trait_mean")
    dat_mean$Genotype <- as.character(dat_mean$Genotype)
    # Prep BLUP data 
    blup_med <- dplyr::left_join(BLUP %>% dplyr::select(c("Genotype", "blup_trait")), dat_median, by = "Genotype") 
    blup_median_mean <- dplyr::left_join(blup_med, dat_mean, by = "Genotype")
    
    p1 <- ggplot2::ggplot(blup_median_mean, ggplot2::aes(x = .data[["blup_trait"]], y = .data[["Trait_median"]])) + 
      ggplot2:: geom_point(na.rm = TRUE) +
      ggplot2::ggtitle(paste0("Trait median vs BLUP: ", t)) + 
      ggplot2::theme_light() +
      ggplot2::stat_smooth(method = "lm", formula = y ~ x, se = FALSE, na.rm = TRUE) +
      ggpubr::stat_cor(method = "pearson", na.rm = TRUE)
    p2 <- ggplot2::ggplot(blup_median_mean, ggplot2::aes(x = .data[["blup_trait"]], y = .data[["Trait_mean"]])) + 
      ggplot2:: geom_point(na.rm = TRUE) +
      ggplot2::ggtitle(paste0("Trait mean vs BLUP: ", t)) + 
      ggplot2::theme_light() +
      ggplot2::stat_smooth(method = "lm", formula = y ~ x, se = FALSE, na.rm = TRUE) +
      ggpubr::stat_cor(method = "pearson", na.rm = TRUE)
    p3 <- ggplot2::ggplot(blup_median_mean, ggplot2::aes(x = .data[["Trait_mean"]], y = .data[["Trait_median"]])) + 
      ggplot2:: geom_point(na.rm = TRUE) +
      ggplot2::ggtitle(paste0("Trait mean vs Trait median: ", t)) + 
      ggplot2::theme_light() +
      ggplot2::stat_smooth(method = "lm", formula = y ~ x, se = FALSE, na.rm = TRUE) +
      ggpubr::stat_cor(method = "pearson", na.rm = TRUE)
    
    g1 <- ggplot2::ggplotGrob(p1)
    g2 <- ggplot2::ggplotGrob(p2)
    g3 <- ggplot2::ggplotGrob(p3)
    g <- cbind(g1, g2, g3, size = "first")
      })
  
  message("Creating correlation plots")
  ml <- gridExtra::marrangeGrob(plot_lst, nrow = 3, ncol = 1)
  ggplot2::ggsave(file.path(stats_plot_dir, "Trait_median_mean_vs_BLUP_corplots.pdf"), ml, width = 15, height = 10)
  }

# Plot uniform
if (is.list(results)) {
  pdf(file = file.path(stats_plot_dir, "DHARMa_model_validation_plots.pdf"), width = 10, height = 6)
  purrr::map(names(results$Tweedie_residual_simulation), function(t) {
    # Extract proportion of zero
    prop_zero <- mean(!results$Outlier_removed_data[[t]], na.rm = TRUE)

    # Make plots 
    p.multi.plot %<a-% {
      graphics::split.screen(c(1,2))
    
      graphics::screen(1)
      DHARMa::plotQQunif(results$Tweedie_residual_simulation[[t]])
      graphics::mtext(paste0("TWEEDIE: ", t," (", round(prop_zero, 3)*100, "% zeros)"), side = 3)
      
      graphics::screen(2)
      DHARMa::plotResiduals(results$Tweedie_residual_simulation[[t]])

    graphics::close.screen(all = TRUE)
  }
  p.multi.plot
  })
graphics.off()
}

#######################################################################
# Save session info
#######################################################################

writeLines(capture.output(sessionInfo()), con = file.path(stats_out_dir, paste("sessionInfo",Sys.Date() ,".txt", sep="")))
message("Done")