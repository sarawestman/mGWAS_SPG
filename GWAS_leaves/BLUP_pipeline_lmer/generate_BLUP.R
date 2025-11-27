#' generate BLUP function
#'
#' @description The goal of generate_BLUP is to run Best Linear Unbiased Predictions.
#' @param dat An input dataset.
#' @param random_effect The random effect column(s). The first column number in vector will be the selected BLUP intercept. Random effect will be included in the model as: (1 | random effect).
#' @param sample_column The sample ID column.
#' @param start_column The start column index for traits.
#' @param fixef_contrasts Sets contrast coding system. If TRUE, compares each fixed effect factor level to a baseline based on the grand mean of all levels (i.e., deviation coding). If FALSE, compares each of the fixed effect factor levels to a reference level based on the first level of a factor (i.e., treatment coding). Notes: If you plan to use the fixed effect estimates using deviation coding, be aware that the level names have been lost and that the estimates require transformation.
#' @param fixed_effect The fixed effect column(s).
#' @param col_convert If TRUE, converts model effect columns to factor type.
#' @keywords BLUP, Pre-GWAS
#' @export

generate_BLUP <- function(dat, random_effect, sample_column, start_column, fixed_effect = NULL, col_convert = TRUE, fixef_contrasts = TRUE) {

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

if (missing(sample_column)) {
  stop("Sample column needs to be provided")
} else if (!is.numeric(sample_column) & sample_column < 1) {
  stop("Sample column needs to be > 0")
}

if (missing(start_column)) {
  stop("Start column needs to be provided")
} else if (!is.numeric(start_column) & start_column < 1) {
  stop("Start column needs to be > 0")
}

if (!is.null(fixed_effect)) {
  if(!is.numeric(fixed_effect) & any(fixed_effect < 1) | length(fixed_effect) == 0) {
    stop("Fixed effect column number(s) must be integers > 0")
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

# Convert the response columns to numeric
dat <- dplyr::mutate(dat, across(seq(start_column, ncol(dat)), as.numeric))

# Use this sanity check to see if all columns (except traits) are consistent throughout the analysis.
bf_start_column <- start_column-1
start_colnames <- colnames(dat)[1:bf_start_column]

#######################################################################
## Outlier removal
#######################################################################
message("Remove outliers")

# Create lmer formula
termlabels <- c(paste0("(1|", colnames(dat)[random_effect], ")"))
if (length(random_effect) > 0 & !is.null(fixed_effect)) { # with fixed effect
  termlabels <- c(colnames(dat)[fixed_effect], termlabels)
}
message(paste("\tlmer model, effect included in formula:", paste(termlabels, collapse = ", ")))

# Find outliers
trait_names <- colnames(dat)[start_column:ncol(dat)]

outliers_residuals <- purrr::map(setNames(trait_names, trait_names), function(t) {
  res <- car::outlierTest(
    lme4::lmer(formula = reformulate(termlabels = termlabels, response = t),
               data = dat, REML = TRUE),
    n.max = nrow(dat),
    cutoff = 0.05
  )
  as.integer(names(res$bonf.p)[res$bonf.p < 0.05])
})

# Store outlier samples in a data frame
outlier_dat <- purrr::map2_dfr(names(outliers_residuals), outliers_residuals, ~ {
  if (length(.y) > 0) {
    dplyr::select(dat, all_of(sample_column), all_of(.x)) %>% dplyr::slice(.y)
  }
})

# Remove outlier samples 
purrr::walk2(names(outliers_residuals), outliers_residuals, ~ {
  if (length(.y) > 0) {
    dat <<- dplyr::mutate(dat, across(all_of(.x), function(column) replace(column, .y, NA)))
  }
})

outlier_removed_dat <- dat

# Sanity check 
if(!identical(start_colnames, colnames(dat)[1:bf_start_column])){
  stop("Column names (before start_column) have changed. The column assumption no longer hold.")
}

#######################################################################
# Normalization 
#######################################################################
message("Check distribution of residuals")
# Run lmer model
lmms_all <- lapply(trait_names, function(t) {
  lme <- lme4::lmer(formula = reformulate(termlabels = termlabels, response = t), data = dat, REML = TRUE)
})
names(lmms_all) <- trait_names

# Make BLUPs
blup <- lapply(trait_names, function(t) {
  lme <- lmms_all[[t]]
  message(paste0("\tBLUP selected (", t, "): ", names(lme4::ranef(lme))[1]))
  blup <- lme4::ranef(lme)[[1]][[1]] 
  data.frame(row.names = lme4::ranef(lme)[[1]] %>% rownames(),
             Genotype = lme4::ranef(lme)[[1]] %>% rownames(),
             blup = blup,
             trait = t)
}) %>% bind_rows %>% tibble::as_tibble() %>% tidyr::spread(trait, blup)

# Remove isSingular and non-varying traits 
message("Check for isSingular and zero variance traits")
trait_isSingular1 <- tibble::tibble(Trait = names(lmms_all), singular = purrr::map_lgl(lmms_all, lme4::isSingular)) %>% dplyr::filter(singular) %>% dplyr::pull(Trait)
traits_zero_var1 <- blup %>% dplyr::summarise_all(var) %>% dplyr::select_if(function(.) . == 0) %>% names()
trait_zerosing1 <- c(traits_zero_var1, trait_isSingular1)

if(length(trait_zerosing1) > 0) {
  lapply(trait_isSingular1, function(t) {
    message(paste0("\tisSingular traits: ", t))
  })
  lapply(traits_zero_var1, function(t) {
    message(paste0("\tZero variance traits: ",t))
  })
  trait_zerosing1 <- unique(trait_zerosing1)
  trait_names <- trait_names[!(trait_names %in% trait_zerosing1)]
  dat <- dat %>% dplyr::select(-all_of(trait_zerosing1))
  lapply(trait_zerosing1, function(t) {
    message(paste0("\tTraits removed: ",t))
  })
  trait_zerosing <- trait_zerosing1
}

# Store shapiro test results
shap_norm <- purrr::map_df(trait_names, function(t) {
  x <- lmms_all[[t]]
  ranef_ls <- lapply(names(lme4::ranef(x)), function(i){ # run shapiro test on each random effect 
    shap_ranef <- shapiro.test(lme4::ranef(x)[[i]][,1])
    df <- tibble(Trait = t,
                 ranef_pval = shap_ranef$p.value,
                 ranef_isNorm = ranef_pval > 0.05)
    df <- setNames(df, c("Trait", paste0(colnames(df)[2:3], "_", i))) # add random effect name to columns  
  }) 
  ranef_df <- Reduce(function(x, y) merge(x , y, by = "Trait"), ranef_ls) # merge random effect list to a data frame  
  shap_resid <- shapiro.test(resid(x))
  resid_df <- tibble::tibble(Trait = t,
                     resid_pval = shap_resid$p.value,
                     resid_isNorm = resid_pval > 0.05)
  full_df <- dplyr::full_join(resid_df, ranef_df, by = "Trait")
})

## Extract isNorm and trait columns 
isNorm_df <- shap_norm %>% dplyr::select("Trait", dplyr::starts_with("ranef_isNorm"), dplyr::starts_with("resid_isNorm")) %>% as.data.frame() %>% tibble::column_to_rownames(var = "Trait")

## Extract trait names that have normally distributed residuals AND random effect residuals
non_trans_trait <- apply(isNorm_df, 1, function(x) sum(x == TRUE) / length(x)) %>% as.data.frame() %>% subset(. == 1) %>% rownames() 
## Extract trait names that do not fulfill the above conditions 
trans_trait <- trait_names[!(trait_names %in% non_trans_trait)]

# Transform those traits
message("Normalize data")
## transform 
dat <- dplyr::mutate(dat, dplyr::across(all_of(trans_trait), ~ {
  output <- purrr::quietly(bestNormalize::orderNorm)(.)
  if (length(output$warnings) > 0) {
    message(paste0("\tWarning: ", dplyr::cur_column(), ": ", output$warnings), appendLF = FALSE)
  }
  output$result$x.t
}))

# Sanity check 
if(!identical(start_colnames, colnames(dat)[1:bf_start_column])){
  stop("Column names (before start_column) have changed. The column assumption no longer hold.")
}

#######################################################################
## Generate BLUP
#######################################################################
message("Calculate BLUPs")
message(paste("\tlmer model, effect included in formula:", paste(termlabels, collapse = ", ")))

# Check all columns for NaN and Inf
dat <- dplyr::mutate(dat, across(seq(start_column, ncol(dat)),
                                 ~ ifelse(is.infinite(.) | is.nan(.), NA, .)))

dat_norm <- dat 

# Fit the model
lmms_all <- lapply(trait_names, function(t) {
  lme <- lme4::lmer(formula = reformulate(termlabels = termlabels, response = t), data = dat, REML = TRUE)
})
names(lmms_all) <- trait_names

blup <- lapply(trait_names, function(t) {
  lme <- lmms_all[[t]]
  message(paste0("\tBLUP selected (", t, "): ", names(lme4::ranef(lme))[1]))
  blup <- lme4::ranef(lme)[[1]][[1]] 
  data.frame(row.names = lme4::ranef(lme)[[1]] %>% rownames(),
             Genotype = lme4::ranef(lme)[[1]] %>% rownames(),
             blup = blup,
             trait = t)
}) %>% bind_rows %>% tibble::as_tibble() %>% tidyr::spread(trait, blup)

# Remove isSingular and non-varying traits 
message("Check for isSingular and zero variance traits")
trait_isSingular2 <- tibble::tibble(Trait = names(lmms_all), singular = purrr::map_lgl(lmms_all, lme4::isSingular)) %>% dplyr::filter(singular) %>% dplyr::pull(Trait)
traits_zero_var2 <- blup %>% dplyr::summarise_all(var) %>% dplyr::select_if(function(.) . == 0) %>% names()
trait_zerosing2 <- c(traits_zero_var2, trait_isSingular2)

if(length(trait_zerosing2) > 0) {
  lapply(trait_isSingular2, function(t) {
    message(paste0("\tisSingular traits: ", t))
  })
  lapply(traits_zero_var2, function(t) {
    message(paste0("\tZero variance traits: ",t))
  })
  trait_zerosing2 <- unique(trait_zerosing2)
  trait_names <- trait_names[!(trait_names %in% trait_zerosing2)]
  lmms_all <- lmms_all[!names(lmms_all) %in% trait_zerosing2]
  blup <- dplyr::select(blup, -all_of(trait_zerosing2))
  lapply(trait_zerosing2, function(t) {
    message(paste0("\tTraits removed: ",t))
  })
  if(exists("trait_zerosing1")){
    trait_zerosing <- c(trait_zerosing1, trait_zerosing2)
  } else {
    trait_zerosing <- trait_zerosing2
  }
}

# Store shapiro test
shap_norm_BLUP <- purrr::map_df(trait_names, function(t) {
  x <- lmms_all[[t]]
  ranef_ls <- lapply(names(lme4::ranef(x)), function(i){ # run shapiro test on each random effect 
    shap_ranef <- shapiro.test(lme4::ranef(x)[[i]][,1])
    df <- tibble(Trait = t,
                 ranef_pval = shap_ranef$p.value,
                 ranef_isNorm = ranef_pval > 0.05)
    df <- setNames(df, c("Trait", paste0(colnames(df)[2:3], "_", i))) # add random effect name to columns  
    }) 
  ranef_df <- Reduce(function(x, y) merge(x , y, by = "Trait"), ranef_ls) # merge random effect list to a data frame  
  shap_resid <- shapiro.test(resid(x))
  resid_df <- tibble::tibble(Trait = t,
                     resid_pval = shap_resid$p.value,
                     resid_isNorm = resid_pval > 0.05)
  full_df <- dplyr::full_join(resid_df, ranef_df, by = "Trait")
  })

# Store results 
if(exists("trait_zerosing")) {
list(
  "Outlier_removed_data" = outlier_removed_dat,
  "Outlier_data" = outlier_dat,
  "Outliers_residuals" = outliers_residuals,
  "Shapiro_raw_data" = shap_norm,
  "BLUP" = blup,
  "lme_data" = lmms_all, 
  "Shapiro_BLUP" = shap_norm_BLUP,
  "Norm_raw_dat" = dat_norm,
  "Transformed_traits" = as.data.frame(trans_trait),
  "Non_transformed_traits" = as.data.frame(non_trans_trait),
  "Removed_traits" = as.data.frame(trait_zerosing)
)
} else{
list(
  "Outlier_removed_data" = outlier_removed_dat,
  "Outlier_data" = outlier_dat,
  "Outliers_residuals" = outliers_residuals,
  "Shapiro_raw_data" = shap_norm,
  "BLUP" = blup,
  "lme_data" = lmms_all, 
  "Shapiro_BLUP" = shap_norm_BLUP,
  "Norm_raw_dat" = dat_norm,
  "Transformed_traits" = as.data.frame(trans_trait),
  "Non_transformed_traits" = as.data.frame(non_trans_trait)
  )
}
}
  