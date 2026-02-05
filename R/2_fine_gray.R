# Title:         Analyses based on Fine-Gray models
# Author:        Eva N. S. Wandall
# Reviewer:      Rune Haubo B. Christensen
# Date:          2026-02-04
#
# Description ------------------------------------------------------------------
# This R-script replicates the central analyses of the diagnostic recurrence
# and pairwise diagnostic shifts.
#
# The analyses are demonstrated using females diagnosed with a personality 
# disorder (ICD-10 code F60) studied for the outcome diagnostic recurrence as an
# example. However, the code can be readily adapted to study different 
# populations or pairwise diagnostic shifts.
#
# Content ----------------------------------------------------------------------
# - Subset and format data
# - Fit Fine-Gray model
# - Testing significance of multi-level categorical covariate (Wald test)
# - Compute population averaged cumulative incidence functions
# - Compute age-specific cumulative incidence function
# - Smoothing of the cumulative incidence curve
# - Compute Wald confidence intervals based on bootstrapping
#
# Requirements -----------------------------------------------------------------
# To run this program, the following data are required:
#
# - data_recurrence: data.table object created in 1_format_data.R.
#
# That is, you should be able to run the following lines (after loading the 
# data.table package):
# data_recurrence[, .(id,
#                    f_first_diagnosis,
#                    f_event_type,
#                    f_sex,
#                    f_age, 
#                    f_calendar_time,
#                    tstop,
#                    i_status)] 
#

# ------------------------------------------------------------------------------
# Packages ---------------------------------------------------------------------
# ------------------------------------------------------------------------------

library(data.table) # For data handling
library(survival)   # For finegray(), coxph(), and survfit()
library(timereg)    # For wald.test()
library(ggplot2)    # For visualizations

# ------------------------------------------------------------------------------
# Subset and format data -------------------------------------------------------
# ------------------------------------------------------------------------------

# Specify the first diagnostic group and later diagnostic group (which will be 
# the same when the outcome is diagnostic recurrence) and sex of the population 
# of interest. 
first_diagnosis <- "F60"
later_diagnosis <- "F60"
sex <- "Female"

# Subset data
analysis_data <- data_recurrence[f_first_diagnosis == first_diagnosis &
                                   f_event_type == later_diagnosis &
                                   f_sex == sex]

# Convert i_status to factor variable
analysis_data[, f_status := factor(i_status)]

# ------------------------------------------------------------------------------
# Fit Fine-Gray model ----------------------------------------------------------
# ------------------------------------------------------------------------------

# Create data for a Fine-Gray model
fg_data <- finegray(Surv(tstop, f_status) ~ ., 
                    data = analysis_data)

# Fit the Fine-Gray model
fg_fit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 
                  f_age + 
                  f_calendar_time, 
                weight =fgwt, 
                data = fg_data, 
                ties = "breslow")

# Display summary of model fit (not reported in the manuscript or supplementary, 
# since sub-distribution hazard ratios should be interpreted with caution as 
# described in eAppendix2.)
summary(fg_fit)

# ------------------------------------------------------------------------------
# Testing significance of multi-level categorical covariate (Wald test) --------
# ------------------------------------------------------------------------------

# Significance test of age assuming four levels (presented in eTable 16)
wald.test(coef = fg_fit$coef,     # Coefficient object from the Fine-Gray model
          vcov = vcov(fg_fit),    # Variance of estimates
          coef.null = c(1, 2, 3)) # Indices set to 0 under the null hypothesis

# ------------------------------------------------------------------------------
# Compute population averaged cumulative incidence functions -------------------
# ------------------------------------------------------------------------------

# Create prediction data set containing all individuals
pred_data_pop <- data.table(f_age = analysis_data$f_age,
                            f_calendar_time = analysis_data$f_calendar_time)

# Fit cumulative incidence functions for each individual 
cif_pop_fit <- survfit(fg_fit, newdata = pred_data_pop)

# Extract predicted cumulative incidences 18 years after first contact
cif_pop_t18 <- summary(cif_pop_fit, times = 18)

# Average over predictions to obtain an estimate of the standardized cumulative 
# incidence at 18 years adjusted to the observed covariate distribution (part 
# of Figure 4)
(std_cif_pop_t18 <- rowMeans(1 - cif_pop_t18$surv))

# ------------------------------------------------------------------------------
# Compute age-specific cumulative incidence function ---------------------------
# ------------------------------------------------------------------------------

# We illustrate the age-specific cumulative incidence estimation using the age
# group 18-24 years. By changing the 'age_level' argument predictions for other
# levels can be obtained. 

# Specify age group of interest 
age_level <- "18-24" 

# Create prediction data set containing all individuals where the covariate of 
# interest (here age) has been set to a fixed level (specified by the age_level
# argument)
pred_data_age <- data.table(f_age = age_level,
                            f_calendar_time = analysis_data$f_calendar_time)

# Fit cumulative incidence functions for each individual 
cif_fit_age <- survfit(fg_fit, newdata = pred_data_age)

# Compute the standardized cumulative incidence function for individuals
# aged 18-24 diagnosed with a personality disorder at first contact by averaging 
# over the predictions
std_cif_age <- rowMeans(1 - cif_fit_age$surv)

# Collect results in a data.table
cif_data_age <- data.table(time = cif_fit_age$time, 
                           cuminc = std_cif_age, 
                           age = age_level)

# ------------------------------------------------------------------------------
# Smoothing of the cumulative incidence curve ----------------------------------
# ------------------------------------------------------------------------------

# The cumulative incidence functions presented in the manuscript are smoothed to 
# comply with the Danish Data Protection Regulations

# Log transformation of time 
cif_data_age[, logtime := log(time)]

# Fit a smoothing spline
cuminc_spline <- cif_data_age[time <= 18, 
                              smooth.spline(x = logtime,
                                            y = cuminc,
                                            spar = 0.85)] # sparsity level

# Extract values from the spline object
cif_data_age[time <= 18, cuminc_smoothed := cuminc_spline$y]

# Plot result with ggplot2 (part of eFigure 7)
ggplot(data = cif_data_age[time <= 18], aes(x = time, y = cuminc_smoothed)) + 
  geom_step() + 
  xlab("Time since first contact (years)") +
  ylab("Cumulative incidence") +
  xlim(0, 18) +
  theme_bw()

# ------------------------------------------------------------------------------
# Compute Wald confidence intervals based on bootstrapping ---------------------
# ------------------------------------------------------------------------------

# We illustrate the procedure for obtaining confidence intervals based on 
# bootstrapping for the age-specific cumulative incidence 18 years after first 
# contact for the individuals aged 18-24 years at first contact. By changing 
# the 'age_level' and 't' arguments in the boot function, predictions for other 
# age levels and time points can be obtained. Alternatively, the boot package 
# (Canty & Ripley, 2022) can be used to compute a bootstrap estimate of the 
# standard error.

# Specification of age level and time point of interest as well as the desired
# number of bootstrap samples:

# - age_level:      Age level of interest:
# - t:              Time point of interest
# - R:              Number of boostrap samples

age_level <- "18-24"
t <- 18
R <- 2

# Function that computed the standardized cumulative incidence at time 't' 
# for individuals in the age group 'age_level' on a bootstrap sample 
boot_std_cuminc_at_t <- function(age_level, t) {
  # Sample data with replacement
  sample_data <-  analysis_data[sample(.N, .N, replace = TRUE)]
  
  # Create data suitable for fitting a Fine-Gray model
  fg_data <- finegray(Surv(tstop, as.factor(f_status)) ~ ., 
                      data = sample_data)
  
  # Fit Fine-Gray model to the bootstrap sample
  fg_fit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 
                    f_age + f_calendar_time, 
                  weight = fgwt, 
                  data = fg_data, 
                  ties = "breslow")
  
  # Number of individuals within each covariate stratum
  covariate_strata <- sample_data[, .N, by = c("f_age", "f_calendar_time")]
  
  # Add identification variable for all stratum
  covariate_strata[, strata_id := factor(1:.N)]
  
  # Modify the age level such that predictions correspond to a fixed level (set
  # by age_level) across all strata
  covariate_strata[, f_age := age_level]
  
  # Predict cumulative incidence at time 't' for all strata
  cif_pred <- survfit(fg_fit, newdata = covariate_strata)
  cif_at_t <- summary(cif_pred, times=t)
  cif_data <- data.table(time = cif_at_t$time, cif_at_t$surv)
  
  # Attach stratum sizes for weighted averaging
  cif_long <- melt(cif_data, 
                   id.vars = "time", 
                   value = "subdist_surv", 
                   variable = "strata_id")
  cif_long[, strata_id := factor(strata_id)]
  cif_long <- merge(cif_long, 
                    covariate_strata[, .(strata_id,N)], 
                    all.x = TRUE, 
                    by = "strata_id")
  
  # Compute standardized estimate at time 't' on bootstrap sample
  std_cif <- cif_long[, sum(N*(1 - subdist_surv))/sum(N)]
  
  # Return:
  std_cif
}

# Seed the random number generator for reproducibility
set.seed(1234)

# Sample data with replacement R times and compute the statistic of interest 
# on each sample
boot_results <- sapply(1:R, function(x) boot_std_cuminc_at_t(age_level, t))

# The estimated cumulative incidence of diagnostic recurrence 18 years after 
# first contact for females aged 18-24 diagnosed with a personality disorder at 
# first contact
(std_cuminc_at_t <- cif_data_age[time < t][.N]$cuminc)

# Bootstrap estimated standard error 
(boot_se <- sd(boot_results))

# 95% Wald confidence interval
(std_cuminc_at_t + c(-1, 1) * qnorm(0.975) * boot_se)

# The estimate and corresponding bootstrap confidence interval are presented in 
# the diagnostic recurrence subsection of the results section.

# ------------------------------------------------------------------------------
# R session info ---------------------------------------------------------------
# ------------------------------------------------------------------------------

sessionInfo()
# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows Server 2022 x64 (build 20348)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8
# [4] LC_NUMERIC=C                    LC_TIME=Danish_Denmark.utf8    
# 
# time zone: Europe/Copenhagen
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] boot_1.3-30       ggplot2_3.5.1     timereg_2.0.5     survival_3.7-0    data.table_1.16.2
# 
# loaded via a namespace (and not attached):
#   [1] vctrs_0.6.5         cli_3.6.3           rlang_1.1.4         generics_0.1.3      labeling_0.4.3     
# [6] glue_1.8.0          colorspace_2.1-1    listenv_0.9.1       future.apply_1.11.2 lava_1.8.0         
# [11] fansi_1.0.6         scales_1.3.0        grid_4.4.1          munsell_0.5.1       tibble_3.2.1       
# [16] lifecycle_1.0.4     numDeriv_2016.8-1.1 compiler_4.4.1      dplyr_1.1.4         codetools_0.2-20   
# [21] pkgconfig_2.0.3     rstudioapi_0.17.0   future_1.33.2       farver_2.1.2        lattice_0.22-6     
# [26] digest_0.6.37       R6_2.5.1            tidyselect_1.2.1    utf8_1.2.4          pillar_1.9.0       
# [31] parallelly_1.45.0   parallel_4.4.1      splines_4.4.1       magrittr_2.0.3      Matrix_1.7-0       
# [36] withr_3.0.0         tools_4.4.1         gtable_0.3.5        globals_0.16.3 
# Program end ------------------------------------------------------------------