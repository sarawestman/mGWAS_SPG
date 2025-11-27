#' generate BLUP function
#'
#' @description The goal of generate_tweedie_BLUP is to run Best Linear Unbiased Predictions. 
#' @param dat An input dataset.
#' @param random_effect The random effect column(s). The first column number in vector will be the selected BLUP intercept. Random effect will be included in the model as: (1 | random effect).
#' @param trait_columns The column index for traits.
#' @param fixef_contrasts Sets contrast coding system. If TRUE, compares each fixed effect factor level to a baseline based on the grand mean of all levels (i.e., deviation coding). If FALSE, compares each of the fixed effect factor levels to a reference level based on the first level of a factor (i.e., treatment coding). Notes: If you plan to use the fixed effect estimates using deviation coding, be aware that the level names have been lost and that the estimates require transformation.
#' @param fixed_effect The fixed effect column(s).
#' @param col_convert If TRUE, converts model effect columns to factor type.
#' @param nsim Specifies the number of simulations in simulateResiduals. If NULL, DHARMa default settings will be used (n = 250). If you want to remove data points, you should at least increase n to > 5 * datasize to decrease the number of "random" outliers, and only get points that are really outside the distribution.
#' @keywords BLUP, Pre-GWAS
#' @export

generate_tweedie_BLUP <- function(dat, random_effect, trait_columns, nsim = NULL, fixed_effect = NULL, col_convert = TRUE, fixef_contrasts = TRUE) {

`%>%` <- dplyr::`%>%`

#######################################################################
## Check and set input arguments 
#######################################################################

# Missing arguments 
if (missing(dat)) {
  stop("Data needs to be provided")
  }

if (missing(random_effect)) {
  stop("Random effect column(s) needs to be provided")
  } else if (!is.numeric(random_effect) & any(random_effect < 1) | length(random_effect) == 0) {
  stop("Random effect needs to be > 0")
    }

if (missing(trait_columns)) {
  stop("Start column needs to be provided")
  } else if (!is.numeric(trait_columns) & any(trait_columns) < 1) {
    stop("Start column needs to be > 0")
    }

if (!is.null(fixed_effect)) {
  if(!is.numeric(fixed_effect) & any(fixed_effect < 1) | length(fixed_effect) == 0) {
    stop("Fixed effect column number(s) must be integers > 0")
  }
  }

if (!is.null(nsim)) {
  if(!is.numeric(nsim) & any(nsim < 1) | length(nsim) == 0) {
    stop("Variable nsim must be a integer > 0")
  }
  }

# Effect column conversion - fixed and random effects 
if (col_convert) {
  message("\tConverting lmer model effect columns to factor type")
  dat <- dplyr::mutate(dat, across(all_of(random_effect), as.factor))
  if (!is.null(fixed_effect)) {
    dat <- dplyr::mutate(dat, across(all_of(fixed_effect), as.factor))
  }
}

# Contrast coding - fixed effects 
if (!is.null(fixed_effect)) {
  if (fixef_contrasts) {
  message("\tContrast coding system: Deviation coding")
    for(i in fixed_effect) {
      if(is.factor(dat[[i]])) {
        contrasts(dat[[i]]) <- contr.sum(dat[[i]] %>% levels() %>% length())
        } else {    
          stop("Contrasts require fixed effects as factor type")
        }
    }
  }
  }

# Extract trait names
trait_names <- colnames(dat)[trait_columns]

# Convert the response columns to numeric
dat <- dplyr::mutate(dat, across(all_of(trait_names), as.numeric))

#######################################################################
## Fit model and simulate residuals 
#######################################################################

# Create formula
termlabels <- c(paste0("(1|", colnames(dat)[random_effect], ")"))
if (length(random_effect) > 0 & !is.null(fixed_effect)) { # with fixed effect
  termlabels <- c(colnames(dat)[fixed_effect], termlabels)
  }
message(paste("Tweedie model, effect included in formula:", paste(termlabels, collapse = ", ")))

# Fit the model  
message("Fit tweedie model")
glmm_tweedie <- purrr::map(setNames(trait_names, trait_names), function(t) {
    fittedModel <- glmmTMB::glmmTMB(reformulate(termlabels = termlabels, response = t), family = glmmTMB::tweedie(link = "log"), data = dat)
    })

# Residual simulations 
message("Simulate residuals")
if (!is.null(nsim)) {
  tweedie_sim <- purrr::map(setNames(trait_names, trait_names), function(t) {
    simulationOutput <- DHARMa::simulateResiduals(fittedModel = glmm_tweedie[[t]], n = nsim, plot = FALSE)
    })
  } else {
  tweedie_sim <- purrr::map(setNames(trait_names, trait_names), function(t) {
    simulationOutput <- DHARMa::simulateResiduals(fittedModel = glmm_tweedie[[t]], n = 250, plot = FALSE)
  })
  }

#######################################################################
## Check for zero-inflated traits 
#######################################################################
message("Check for zero-inflated traits")

# Check for zero-inflation
tweedie_zeroinf <- purrr::map_df(names(tweedie_sim), function(t) {
  res <- DHARMa::testZeroInflation(tweedie_sim[[t]], plot = FALSE)
  tibble(Trait = t, 
         ratioObsSim = res$statistic[[1]],
         ZeroInf_pval = res$p.value)
  })

traits_noInfl <- tweedie_zeroinf %>% subset(ZeroInf_pval > 0.05) %>% dplyr::pull(Trait)
traits_Infl <- tweedie_zeroinf$Trait[!tweedie_zeroinf$Trait %in% traits_noInfl]

if (length(traits_Infl) > 0) {
  lapply(traits_Infl, function(t) {
    message(paste0("\tZero-inflated tweedie trait: ", t))
    })
  }

#######################################################################
## Outlier removal
#######################################################################
message("Remove outliers")

# Find outliers
outliers_residuals <- purrr::map(setNames(trait_names, trait_names), function(t) {
  res <- DHARMa::outliers(tweedie_sim[[t]], lowerQuantile = 0, upperQuantile = 1) # DHARMA default settings
  })

# Store outlier samples in a data frame
outlier_dat <- purrr::map2_dfr(names(outliers_residuals), outliers_residuals, ~ {
  if (length(.y) > 0) {
    dplyr::select(dat, all_of(names(tweedie_sim)), all_of(.x)) %>% dplyr::slice(.y)
    }
  })

# Remove outlier samples 
purrr::walk2(names(outliers_residuals), outliers_residuals, ~ {
  if (length(.y) > 0) {
    dat <<- dplyr::mutate(dat, across(all_of(.x), function(column) replace(column, .y, NA)))
    }
  })

outlier_removed_dat <- dat

#######################################################################
# Validate model
#######################################################################
message("Model assumption stats")

# Check all columns for NaN and Inf
dat <- dplyr::mutate(dat, across(all_of(trait_names),
                                 ~ ifelse(is.infinite(.) | is.nan(.), NA, .)))

# Fit the model
glmm_tweedie <- purrr::map(setNames(trait_names, trait_names), function(t) {
    fittedModel <- glmmTMB::glmmTMB(reformulate(termlabels = termlabels, response = t), family = glmmTMB::tweedie(link = "log"), data = dat)
  })

# Make residual simulations 
if (!is.null(nsim)) {
  tweedie_sim <- purrr::map(setNames(trait_names, trait_names), function(t) {
    simulationOutput <- DHARMa::simulateResiduals(fittedModel = glmm_tweedie[[t]], n = nsim, plot = FALSE)
    })
  } else {
    tweedie_sim <- purrr::map(setNames(trait_names, trait_names), function(t) {
      simulationOutput <- DHARMa::simulateResiduals(fittedModel = glmm_tweedie[[t]], n = 250, plot = FALSE)
      })
  }

# Model assumption stats 
testRes_stats <- purrr::map_df(trait_names, function(t) {
  invisible(capture.output(res <- DHARMa::testResiduals(tweedie_sim[[t]], plot = FALSE)))
  tibble(Trait = t,
         Uni_tweedie = res$uniformity$p.value,
         Disp_tweedie = res$dispersion$p.value,
         Outl_tweedie = res$outlier$p.value,)
  })

#######################################################################
## Generate BLUP
#######################################################################
message("Calculate BLUPs")

# Make BLUPs 
blup <- lapply(trait_names, function(t) {
  my_tweedie <- glmm_tweedie[[t]]
  message(paste0("\tBLUP selected (", t, "): ", names(glmmTMB::ranef(my_tweedie)[[1]])))
  blup <- glmmTMB::ranef(my_tweedie)[[1]][[1]][,1]
  data.frame(row.names = glmmTMB::ranef(my_tweedie)[[1]][[1]] %>% rownames(),
             Genotype = glmmTMB::ranef(my_tweedie)[[1]][[1]] %>% rownames(),
             blup = blup,
             trait = t)
  }) %>% bind_rows %>% tibble::as_tibble() %>% tidyr::spread(trait, blup)

# Remove non-varying traits 
message("Check for zero variance traits")
traits_zero_var <- blup %>% dplyr::summarise_all(var) %>% dplyr::select_if(function(.) . == 0) %>% names()

if (length(traits_zero_var) > 0) {
  lapply(traits_zero_var, function(t) {
    message(paste0("\tZero variance traits: ",t))
  })
  trait_names <- trait_names[!(trait_names %in% traits_zero_var)]
  glmm_tweedie <- glmm_tweedie[!names(glmm_tweedie) %in% traits_zero_var]
  blup <- dplyr::select(blup, -all_of(traits_zero_var))
  lapply(traits_zero_var, function(t) {
    message(paste0("\tTraits removed: ",t))
  })
}

# Store shapiro test
shap_BLUP <- purrr::map_df(trait_names, function(t) {
  shap_blup <- shapiro.test(blup[[t]])
  tibble(Trait = t,
         BLUP_pval = shap_blup$p.value,
         BLUP_isNorm = BLUP_pval > 0.05)
  })

## Extract isNorm and trait columns 
not_norm_traits <- shap_BLUP %>% subset(BLUP_isNorm == FALSE) %>% dplyr::pull(Trait)
norm_traits <- trait_names[!(trait_names %in% not_norm_traits)]

# Store results 
if(exists("traits_zero_var")) {
  list(
    "Outlier_removed_data" = outlier_removed_dat,
    "Outlier_data" = outlier_dat,
    "Outliers_residuals" = outliers_residuals,
    "testResidual_stats" = testRes_stats,
    "BLUP" = blup,
    "Shapiro_BLUP" = shap_BLUP,
    "Tweedie_residual_simulation" = tweedie_sim,
    "glmm_tweedie" = glmm_tweedie,
    "Zero_inflated_traits" = as.data.frame(traits_Infl),
    "Not_zero_inflated_traits" = as.data.frame(traits_noInfl),
    "Removed_traits" = as.data.frame(traits_zero_var),
    "NotNormDist_BLUP_traits" = as.data.frame(not_norm_traits),
    "NormDist_BLUP_traits" = as.data.frame(norm_traits)
    )
  } else {
    list(
      "Outlier_removed_data" = outlier_removed_dat,
      "Outlier_data" = outlier_dat,
      "Outliers_residuals" = outliers_residuals,
      "testResidual_stats" = testRes_stats,
      "BLUP" = blup,
      "Shapiro_BLUP" = shap_BLUP,
      "Tweedie_residual_simulation" = tweedie_sim,
      "glmm_tweedie" = glmm_tweedie,
      "Zero_inflated_traits" = as.data.frame(traits_Infl),
      "Not_zero_inflated_traits" = as.data.frame(traits_noInfl),
      "NotNormDist_BLUP_traits" = as.data.frame(not_norm_traits),
      "NormDist_BLUP_traits" = as.data.frame(norm_traits)
    )
  }
}
