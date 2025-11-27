#' generate Qst function
#'
#' @description The goal of generate_Qst is to estimate trait Qst (as described in Mähler et al.(2017)) and the confidence interval (using bootstrap).
#' @description Mähler N, Wang J, Terebieniec BK, Ingvarsson PK, Street NR, Hvidsten TR (2017) Gene co-expression network connectivity is an important determinant of selective constraint. PLoS Genet 13(4): e1006402.
#' @param dat An input dataset. 
#' @param random_effect The random effect column(s). The first column number in vector will be the selected BLUP intercept. Random effect will be included in the model as: (1 | random effect).
#' @param start_column The start column index for traits.
#' @param fixef_contrasts Sets contrast coding system. If TRUE, compares each fixed effect factor level to a baseline based on the grand mean of all levels (i.e., deviation coding). If FALSE, compares each of the fixed effect factor levels to a reference level based on the first level of a factor (i.e., treatment coding). Notes: If you plan to use the fixed effect estimates using deviation coding, be aware that the level names have been lost and that the estimates require transformation.
#' @param fixed_effect The fixed effect column(s).
#' @param col_convert If TRUE, converts model effect columns to factor type.
#' @param bootstr_R Size of the posterior sample from the Bayesian bootstrap (default: 1000).
#' @keywords Qst, bootstrap
#' @export

generate_Qst <- function(dat, random_effect, start_column, fixed_effect = NULL, col_convert = TRUE, fixef_contrasts = TRUE, bootstr_R = 1000) {

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
## Check model assumptions 
#######################################################################
message("Check lmer model assumptions")

# Create lmer formula
termlabels <- c(paste0("(1|", colnames(dat)[random_effect], ")"))
if (length(random_effect) > 0 & !is.null(fixed_effect)) { # with fixed effect
  termlabels <- c(colnames(dat)[fixed_effect], termlabels)
}
message(paste("\tlmer model, effect included in formula:", paste(termlabels, collapse = ", ")))

trait_names <- colnames(dat)[start_column:ncol(dat)]

# Check all columns for NaN and Inf
dat <- dplyr::mutate(dat, across(seq(start_column, ncol(dat)),
                                 ~ ifelse(is.infinite(.) | is.nan(.), NA, .)))
# Fit the model
lmms_all <- lapply(trait_names, function(t) {
  lme <- lme4::lmer(formula = reformulate(termlabels = termlabels, response = t), data = dat, REML = TRUE)
  })
names(lmms_all) <- trait_names

# Calculate BLUPs
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
trait_isSingular <- tibble::tibble(Trait = names(lmms_all), singular = purrr::map_lgl(lmms_all, lme4::isSingular)) %>% dplyr::filter(singular) %>% dplyr::pull(Trait)
traits_zero_var <- blup %>% dplyr::summarise_all(var) %>% dplyr::select_if(function(.) . == 0) %>% names()
trait_zerosing <- c(traits_zero_var, trait_isSingular)

if(length(trait_zerosing) > 0) {
  lapply(trait_isSingular, function(t) {
    message(paste0("\tisSingular traits: ", t))
  })
  lapply(traits_zero_var, function(t) {
    message(paste0("\tZero variance traits: ",t))
  })
  trait_zerosing <- unique(trait_zerosing)
  trait_names <- trait_names[!(trait_names %in% trait_zerosing)]
  lmms_all <- lmms_all[!names(lmms_all) %in% trait_zerosing]
  blup <- dplyr::select(blup, -all_of(trait_zerosing))
  lapply(trait_zerosing, function(t) {
    message(paste0("\tTraits removed: ",t))
  })
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

# Sanity check 
if(!identical(start_colnames, colnames(dat)[1:bf_start_column])){
  stop("Column names (before start_column) have changed. The column assumption no longer hold.")
}
#######################################################################
## Generate Qst
#######################################################################
message("Calculate Qst")
message(paste("\tBootstrap posterior sample size:", bootstr_R))

# Qst function 
EstimateQST <- function(dat) {
  qst.calc <- lme4::lmer(formula = reformulate(termlabels = termlabels, response = "trait"), data = dat, REML = TRUE)
  temp <- as.data.frame(lme4::VarCorr(qst.calc))
  pops.var <- temp[2, 5]
  clones.var <- temp[1, 5]
  Qst <- pops.var/(pops.var + 2*clones.var)
  return(Qst)
  }

# Generate Qst table 
qst_table <- purrr::map_df(trait_names, function(t) {
  message("\t","\t",t)
  dat$trait <- dat[[t]]
  bayes.output <- bayesboot::bayesboot(dat, EstimateQST, R = bootstr_R)
  bayes.summary <- unlist(summary(bayes.output))
  tibble::tibble(Trait = t,
                 Qst = bayes.summary[19],
                 Qst_Lower_95percent_CI = bayes.summary[23],
                 Qst_Upper_95percent_CI = bayes.summary[27])
  })

# Store results 
if(exists("trait_zerosing")) {
  list(
    "BLUP" = blup,
    "lme_data" = lmms_all, 
    "Shapiro_BLUP" = shap_norm_BLUP,
    "Qst_table" = qst_table,
    "Removed_traits" = as.data.frame(trait_zerosing)
  )
} else{
  list(
    "BLUP" = blup,
    "lme_data" = lmms_all, 
    "Shapiro_BLUP" = shap_norm_BLUP,
    "Qst_table" = qst_table
    )
  }
}
