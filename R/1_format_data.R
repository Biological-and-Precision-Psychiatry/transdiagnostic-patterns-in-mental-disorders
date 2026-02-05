# Title:         Create data for multi-state survival models
# Author:        Eva N. S. Wandall
# Reviewer:      Rune Haubo B. Christensen
# Date:          2026-02-04
#
#
# Description ------------------------------------------------------------------
# In this R-script we create the data used to analyse four aspects of 
# disease development:
# - The number of psychiatric contacts
# - The number of diagnostic shifts
# - The cumulative incidence of diagnostic recurrence
# - The cumulative incidence of pairwise diagnostic shifts
#
# Content ----------------------------------------------------------------------
# - Construct data for the number of contacts
# - Construct data for the number of diagnostic shifts
# - Construct data for diagnostic recurrence
# - Construct data for pairwise diagnostic shifts
#
# Requirements -----------------------------------------------------------------
# To run this program, input data with the following structure are required:
#
# - data: Individual level data (one record per person per psychiatric contact)
#   stored as a data.table object.
#   - id: Identification number of the individual
#   - d_in: Date of first psychiatric contact
#   - d_exit: Date of emigration, administrative censoring or death
#   - i_death: Binary indicator variable denoting whether a death occurred at 
#     time d_exit. Takes the value 1 if a death is observed, and 0 otherwise.
#   - f_first_diagnosis: Diagnostic group representing the first contact
#     Factor with levels: F10-19, F20, F21-29, F30-31, F32-33, F40-42, F60
#   - f_sex: Sex of the individual. 
#     Factor with levels: Male, Female
#   - f_age: Age group of the individual at first contact.
#     Factor with levels: 18-24, 25-39, 40-59, 60+ (measured in years)
#     If fewer than 10 events are observed in one of the levels, the number of 
#     levels is reduced, cf. eAppendix 2 and eTable 5.
#   - f_calendar_time: Calendar time period at first contact for the individual
#     Factor with levels: 2000-2004, 2005-2009, 2010-2014, 2015-2019
#     If fewer than 10 events are observed in one of the levels, the number of 
#     levels is reduced, cf. eAppendix 2 and eTable 6.
#   - d_start: Start date of current psychiatric contact
#   - f_diagnosis: Diagnostic group representing the current contact
#     Factor with levels: F10-19, F20, F21-29, F30-31, F32-33, F40-42, F60
# 
#
# That is, you should be able to run the following lines (after loading the 
# data.table package):
# data[, .(id,
#          d_in,
#          d_exit,
#          i_death,
#          f_first_diagnosis, 
#          f_sex, 
#          f_age, 
#          f_calendar_time,
#          d_start, 
#          f_diagnosis)]
# 

# ------------------------------------------------------------------------------
# Packages ---------------------------------------------------------------------
# ------------------------------------------------------------------------------

library(data.table) # For data handling

# ------------------------------------------------------------------------------
# Construct data for the number of contacts ------------------------------------
# ------------------------------------------------------------------------------

# Initiate data
data_n_contacts <- copy(data)

# Sort data by individual and start date of contacts
setorder(data_n_contacts, "id", "d_start")

# For each individual, define the end of a contact using the start date of the
# next observed contact
data_n_contacts[, d_stop := shift(d_start, 1, NA, "lead"), by = id]

# Assign status at the end of each contact
# - i_status = 0: censoring
# - i_status = 1: new contact
# - i_status = 2: death
data_n_contacts[!is.na(d_stop), i_status := 1]
data_n_contacts[is.na(d_stop), i_status := 2*i_death]

# Set stop time to the censoring or death date for terminal intervals
data_n_contacts[is.na(d_stop), d_stop := d_exit]

# Convert start and stop dates to numeric (years since first contact)
data_n_contacts[, tstart := as.numeric(as.Date(d_start) - as.Date(d_in))/365.25]
data_n_contacts[, tstop := as.numeric(as.Date(d_stop) - as.Date(d_in))/365.25]

# Add 0.5 day to all censoring or death stop dates to avoid introducing NA's 
# in the model. This does not change the estimation of the number of contacts but 
# it is required for the unbiased estimation of the competing event of death.
data_n_contacts[, tstop := fifelse(i_status %in% c(0, 2), tstop + 0.5/365.25, tstop)]

# ------------------------------------------------------------------------------
# Construct data for the number of diagnostic shifts ---------------------------
# ------------------------------------------------------------------------------

# Initiate data
data_n_shifts <- copy(data)

# Sort data by individual and start date of contacts
setorder(data_n_shifts, "id", "d_start")

# Create a run identifier that increments at each diagnostic change within an
# individual
data_n_shifts[, n_shift := rleid(f_diagnosis) - 1, by = "id"]

# Collapse each run to its first observation
data_n_shifts <- data_n_shifts[, .SD[1], by =c("id", "n_shift")]

# For each individual, define end of diagnostic group using the next observed 
# start date
data_n_shifts[, d_stop := shift(d_start, 1, NA, "lead"), by = id] 

# Assign status at the end of each contact
# - i_status = 0: censoring
# - i_status = 1: diagnostic shift
# - i_status = 2: death
data_n_shifts[!is.na(d_stop), i_status := 1]
data_n_shifts[is.na(d_stop), i_status := 2*i_death]

# Set stop time to the censoring or death date for terminal intervals
data_n_shifts[is.na(d_stop), d_stop := d_exit]

# Convert start and stop dates to numeric (years since first contact)
data_n_shifts[, tstart := as.numeric(as.Date(d_start) - as.Date(d_in))/365.25]
data_n_shifts[, tstop := as.numeric(as.Date(d_stop) - as.Date(d_in))/365.25]

# Add 0.5 day to all censoring or death stop dates to avoid introducing NA's 
# in the model. This does not change the estimation of the number of diagnostic 
# shifts but it is required for the unbiased estimation of the competing event 
# of death.
data_n_shifts[, tstop := fifelse(i_status %in% c(0,2), tstop + 0.5/365.25, tstop)]

# ------------------------------------------------------------------------------
# Construct data for diagnostic recurrence -------------------------------------
# ------------------------------------------------------------------------------

# Restrict to contacts where the diagnostic group of a contact is the same as the 
# diagnostic group recorded at first contact
data_recurrence <- data[f_diagnosis == f_first_diagnosis]

# Sort data by individual and start date of contacts
setorder(data_recurrence, "id", "d_start")

# Add date of next contact and keep only the first row for each individual
data_recurrence[, d_stop := shift(d_start, 1, NA, "lead"), by = id]
data_recurrence <- data_recurrence[, .SD[1], by = id]

# Assign status at the end of each contact
# - i_status = 0: censoring
# - i_status = 1: diagnostic recurrence
# - i_status = 2: death
data_recurrence[!is.na(d_stop), i_status := 1]
data_recurrence[is.na(d_stop), i_status := 2*i_death]

# Set stop date to d_exit for individuals with no observed diagnostic recurrence
data_recurrence[is.na(d_stop), d_stop := d_exit]

# Convert start and stop dates to numeric (years since first contact)
data_recurrence[, tstart := as.numeric(as.Date(d_start) - as.Date(d_in))/365.25]
data_recurrence[, tstop := as.numeric(as.Date(d_stop) - as.Date(d_in))/365.25]

# Add 0.5 day to all censoring or death stop dates to avoid introducing NA's 
# in model. This does not change the estimate of the cumulative incidence of 
# diagnostic recurrence but it is required for the unbiased estimation of the 
# competing event of death.
data_recurrence[, tstop := fifelse(i_status %in% c(0, 2), tstop + 0.5/365.25, tstop)]

# Add endpoint identifier for the event being analyzed
data_recurrence[, f_event_type := factor(f_first_diagnosis)]

# ------------------------------------------------------------------------------
# Construct data for pairwise diagnostic shifts --------------------------------
# ------------------------------------------------------------------------------

# To run this code one must first specify the first and later diagnostic group 
# Here, we demonstrate using the shift from F20 to F32-33

# Specify first diagnostic group
diag1 <- "F20"

# Specify later diagnostic group
diag2 <- "F32-33"

# Restrict to contacts where diag1 and diag2 are given
data_pairwise_shift <- data[f_first_diagnosis == diag1 & 
                              f_diagnosis %in% c(diag1, diag2)]

# Sort data by individual and start date of contacts
setorder(data_pairwise_shift, cols = "id", "d_start")

# Select only first contact and possibly first contact with the later diagnostic 
# group for each individual
data_pairwise_shift <- data_pairwise_shift[, .SD[1], by = c("id", "f_diagnosis")]

# Add date of diagnostic shift and keep only the first row for each 
# individual
data_pairwise_shift[, d_stop := shift(d_start, 1, NA, "lead"), by = id]
data_pairwise_shift <- data_pairwise_shift[, .SD[1], by = id]

# Assign status at the end of each contact
# - i_status = 0: censoring
# - i_status = 1: diagnostic shift
# - i_status = 2: death
data_pairwise_shift[!is.na(d_stop), i_status := 1]
data_pairwise_shift[is.na(d_stop), i_status := 2*i_death]

# Set stop date to d_exit for individuals with no observed pairwise diagnostic 
# shift
data_pairwise_shift[is.na(d_stop), d_stop := d_exit]

# Convert start and stop dates to numeric (years since first contact)
data_pairwise_shift[, tstart := as.numeric(as.Date(d_start) - as.Date(d_in))/365.25]
data_pairwise_shift[, tstop := as.numeric(as.Date(d_stop) - as.Date(d_in))/365.25]

# Add 0.5 day to all censoring or death stop dates to avoid introducing NA's 
# in model. This does not change the estimate of the cumulative incidence of
# the pairwise diagnostic shift but it is required for the unbiased estimation of 
# the competing event of death.
data_pairwise_shift[, tstop := fifelse(i_status %in% c(0, 2), tstop + 0.5/365.25, tstop)]

# Add endpoint identifier for the event being analyzed
data_pairwise_shift[, f_event_type := factor(diag2)]

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
#   [1] data.table_1.16.2
# 
# loaded via a namespace (and not attached):
#   [1] compiler_4.4.1    tools_4.4.1       rstudioapi_0.17.0
# Program end ------------------------------------------------------------------