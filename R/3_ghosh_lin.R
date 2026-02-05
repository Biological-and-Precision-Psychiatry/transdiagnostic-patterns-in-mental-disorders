# Title:         Analyses based on Fine-Gray models
# Author:        Eva N. S. Wandall
# Reviewer:      Rune Haubo B. Christensen
# Date:          2026-02-04
#
# Description ------------------------------------------------------------------
# This R-script replicates the central analyses of the number of psychiatric
# contacts and number of diagnostic shifts.
#
# The analyses are demonstrated using individuals diagnosed with an other 
# psychotic disorder than schizophrenia (ICD-10 codes F21-29) studied for the 
# outcome number of diagnostic shift. However, the code can be readily adapted 
# to study different populations or number of psychiatric contacts.
#
# Content ----------------------------------------------------------------------
# - Subset data
# - Break ties
# - Fit Ghosh-Lin model
# - Extract mean ratios
# - Testing significance of multi-level categorical covariate (Wald test)
# - Compute age-specific number of events
# - Smoothing of the mean cumulative function of the number of events
# - Compute Wald confidence intervals based on bootstrapping
#
# Requirements -----------------------------------------------------------------
# To run this program, the following data are required:
#
# - data_n_shifts: data.table object created in 1_format_data.R.
#
# If the outcome is the number of contacts, data_n_contacts (also created in 
# 1_format_data.R) is required.
#
# That is, you should be able to run the following line (after loading the 
# data.table package):
# data_n_shifts[, .(id,
#                    f_first_diagnosis,
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
library(ggplot2)    # For visualizations
library(mets)       # For recreg()
library(timereg)    # For tie.breaker() and wald.test()

# ------------------------------------------------------------------------------
# Subset data ------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Specify diagnostic group of interest
first_diagnosis <- "F21-29"

# Restrict data to individuals with the first diagnostic group of interest
analysis_data <- data_n_shifts[f_first_diagnosis == first_diagnosis]

# ------------------------------------------------------------------------------
# Break ties -------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Seed the random number generator for reproducibility
set.seed(1234)

# Mark all events (both events of interest and competing events)
analysis_data[, i_event := fifelse(i_status != 0, 1, 0)]

# Break ties with tie.breaker from package timereg
data_no_ties <- tie.breaker(analysis_data, 
                            start = "tstart", 
                            stop = "tstop", 
                            status = "i_event", 
                            id = "id")

# ------------------------------------------------------------------------------
# Fit Ghosh-Lin model ----------------------------------------------------------
# ------------------------------------------------------------------------------

# Fit a Ghosh-Lin model
gl_fit <- recreg(Event(tstart, tstop, i_status) ~ 
                   f_age + 
                   f_sex + 
                   f_calendar_time + 
                   cluster(id), 
                 data = data_no_ties, 
                 cause = 1, 
                 cens.code = 0, 
                 death.code = 2)

# Summary of model
summary(gl_fit)

# ------------------------------------------------------------------------------
# Extract mean ratios ----------------------------------------------------------
# ------------------------------------------------------------------------------

# Mean ratios (presented in eTable 3)
exp(coef(gl_fit))

# ------------------------------------------------------------------------------
# Testing significance of multi-level categorical covariate (Wald test) --------
# ------------------------------------------------------------------------------

# Significance test of age assuming four levels (presented in eTable 3)
wald.test(coef = gl_fit$coef,     # Coefficient object from the Ghosh-Lin model
          Sigma = vcov(gl_fit),   # Variance of estimates
          coef.null = c(1, 2, 3)) # Indices set to 0 under the null hypothesis

# ------------------------------------------------------------------------------
# Compute age-specific number of events ----------------------------------------
# ------------------------------------------------------------------------------

# We illustrate the age-specific number of events using the age group 18-24 
# years. By changing the 'age_level' argument predictions for other covariate 
# levels can be obtained. Sex-specific predictions can also be  obtained by 
# changing which covariate is fixed in the 'covariate_strata' data.table.

age_level <- "18-24"

# Count number of individuals within each covariate stratum
strata_variables <- c("f_sex", "f_age", "f_calendar_time")
covariate_strata <- analysis_data[, .SD[1], by = "id"][, .N, by = strata_variables]

# Create a unique identifier for each stratum
covariate_strata[, id := .I]

# Modify the covariate of interest (here age) such that predictions correspond 
# to a fixed level (set by age_level) across all strata
covariate_strata[, f_age := factor(age_level, 
                                   levels =levels(analysis_data$f_age))]

# Obtain predictions of the mean cumulative event function for each stratum
pred <- predict(gl_fit, covariate_strata)

# Convert predictions to data.table
# - first column corresponds to time points
# - remaining columns corresponds to predictions for covariate strata
pred_data <- data.table(time = pred$time, 
                        as.data.table(t(pred$cumhaz)))

# Reshape predictions to long format
pred_long <- melt(pred_data, 
                  id.vars = "time", 
                  value = "cumhaz")

# Extract stratum id from prediction column names
pred_long[, id := as.integer(sub("V","",variable))]

# Attach stratum sizes for weighted averaging
pred_long <- merge(pred_long, 
                   covariate_strata[, .(id, N)], 
                   all.x= TRUE, 
                   by = "id")

# Compute the standardized cumulative function of the number of diagnostic shifts 
# for individuals aged 18-24 diagnosed with an other psychotic disorder than 
# schizophrenia at first contact by weighting the predictions
std_mcf <- pred_long[, .(mcf = sum(N * cumhaz)/sum(N)), by = time]

# ------------------------------------------------------------------------------
# Smoothing of the mean cumulative function of the number of events ------------
# ------------------------------------------------------------------------------

# The mean cumulative function of the number of events are smoothed to comply 
# with the Danish Data Protection Regulations.

# Log transformation of time 
std_mcf[, logtime := log(time)]

# Fit a smoothing spline
# The sparsity level for the number of contacts was 0.85
mcf_spline <- std_mcf[time <= 18,
                      smooth.spline(x = logtime,
                                    y = mcf,
                                    spar = 0.7)] # sparsity level

# Extract values from the spline object
std_mcf[time <= 18, mcf_smoothed := predict(mcf_spline, logtime)$y]

# Plot result with ggplot2 (part of Figure 3)
ggplot(std_mcf[time <= 18], aes(x = time, y = mcf_smoothed)) + 
  geom_step() + 
  xlab("Time since first contact (years)") +
  ylab("Expected number of diagnostic shifts") +
  xlim(0, 18) +
  theme_bw()

# ------------------------------------------------------------------------------
# Wald confidence intervals based on bootstrapping -----------------------------
# ------------------------------------------------------------------------------

# We demonstrate the procedure for obtaining confidence intervals based on 
# bootstrapping for the age-specific number of cumulative events at 18 years 
# with the age group 18-24 years. By changing the 'age_level' and 't' argument 
# predictions for other covariate levels or time point can be obtained. 
# Sex-specific predictions can also be obtained by changing which covariate is 
# fixed in the 'covariate_strata' data.table in the boot_std_mcf_at_t function.

# Specification of covariate level and time point of interest as well as the 
# desired number of bootstrap samples:

# - age_level:      Age level of interest:
# - t:              Time point of interest
# - R:              Number of boostrap samples

age_level <- "18-24" 
t <- 18
R <- 2

# Function that compute the standardized mean number of events at time 't' for 
# individuals of the 'age_level' on a bootstrap sample of data
boot_std_mcf_at_t <- function(age_level, time) {
  # Sample data with replacement
  ids <- sample(analysis_data[, .SD[1], id]$id, replace = TRUE)
  sample_data <- rbindlist(
    lapply(seq_along(ids), function(x) {
      analysis_data[id == ids[x]][, new_id := x]
    }))
  
  # Break ties in the bootstrap sample
  sample_data[, i_event := ifelse(i_status != 0, 1,0)]
  sample_no_ties <- tie.breaker(sample_data, 
                                start = "tstart", 
                                stop = "tstop",
                                status = "i_event", 
                                id = "new_id")
  
  # Fit the Gosh-Lin model to the bootstrap sample
  gl_fit <- recreg(Event(tstart, tstop, i_status) ~ 
                     f_sex + 
                     f_age + 
                     f_calendar_time + 
                     cluster(new_id),  
                   data = sample_no_ties, 
                   cause = 1, 
                   cens.code = 0, 
                   death.code = 2)
  
  # Number of individuals within each covariate stratum
  covariate_strata <- sample_data[, .SD[1], by = "new_id"][
    , .N, by = strata_variables]
  
  # Add identification variable for all stratum
  covariate_strata[, new_id := 1:.N]
  
  # Modify the age level such that predictions correspond to a fixed level (set
  # by age_level) across all strata
  covariate_strata[, f_age := factor(age_level, 
                                     levels =levels(analysis_data$f_age))]
  
  # Predict expected number of events for all strata
  pred <- predict(gl_fit, covariate_strata)
  pred_data <- data.table(time = pred$time, as.data.table(t(pred$cumhaz)))
  
  # Attach stratum sizes for weighted averaging
  pred_long <- melt(pred_data, id.vars = "time", value = "cumhaz")
  pred_long[, new_id := as.integer(substr(variable, 2,3))]
  pred_long <- merge(pred_long, 
                     covariate_strata[, .(new_id,N)], 
                     all.x = TRUE, 
                     by = "new_id")
  
  # Compute standardized estimate at time 't' in bootstrap sample
  std_mcf <- pred_long[time < t, sum(N*cumhaz)/sum(N), by = time][.N][[2]]
  
  # Return
  std_mcf
}

# Sample data with replacement R times and compute the statistic of interest 
# on each sample of data
boot_results <- sapply(1:R, function(x) boot_std_mcf_at_t(age_level, t))

# Expected number of diagnostic shifts for individuals aged 18-24 diagnosed with 
# an other psychotic disorder than schizophrenia at first contact 18 years after
# first contact
(std_mcf_at_t <- std_mcf[time < t][.N]$mcf)

# Bootstrap estimated standard error 
(boot_se <- sd(boot_results))

# 95% Wald confidence interval
(std_mcf_at_t + c(-1, 1) * qnorm(0.975) * boot_se)

# The estimate and corresponding bootstrap confidence intervals are reported in 
# the number of diagnostic shifts subsection of the results section as well as 
# in eTable 4 and Figure 3.

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
#   [1] mets_1.3.4        timereg_2.0.5     survival_3.7-0    ggplot2_3.5.1     data.table_1.16.2
# 
# loaded via a namespace (and not attached):
#   [1] Matrix_1.7-0        gtable_0.3.5        future.apply_1.11.2 dplyr_1.1.4         compiler_4.4.1     
# [6] Rcpp_1.0.13         tidyselect_1.2.1    parallel_4.4.1      globals_0.16.3      splines_4.4.1      
# [11] scales_1.3.0        lattice_0.22-6      R6_2.5.1            labeling_0.4.3      generics_0.1.3     
# [16] future_1.33.2       tibble_3.2.1        munsell_0.5.1       pillar_1.9.0        rlang_1.1.4        
# [21] utf8_1.2.4          cli_3.6.3           withr_3.0.0         magrittr_2.0.3      digest_0.6.37      
# [26] grid_4.4.1          mvtnorm_1.2-5       rstudioapi_0.17.0   lifecycle_1.0.4     lava_1.8.0         
# [31] vctrs_0.6.5         glue_1.8.0          farver_2.1.2        numDeriv_2016.8-1.1 listenv_0.9.1      
# [36] codetools_0.2-20    parallelly_1.45.0   fansi_1.0.6         colorspace_2.1-1    tools_4.4.1        
# [41] pkgconfig_2.0.3 
# Program end ------------------------------------------------------------------