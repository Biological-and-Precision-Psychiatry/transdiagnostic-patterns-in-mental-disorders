

# Nationwide transdiagnostic age and sex patterns in mental disorder disease developments


- Authors: [Eva N. S. Wandall](https://github.com/evaninawandall), 
[Per K. Andersen](https://researchprofiles.ku.dk/da/persons/per-kragh-andersen/), 
[Rune H. B. Christensen](https://github.com/runehaubo), 
[Michael E. Benros](https://researchprofiles.ku.dk/da/persons/jrd819-jrd819/)
- A scientific paper under review in [JAMA Psychiatry](https://jamanetwork.com/journals/jamapsychiatry)


This repository contains R-code to replicate the results presented in the paper.


## R-scripts

The R-scripts listed below are available in the `R` folder in this repository. The files contain code that reproduces key analyses of the paper. The data to run the code cannot be made publicly available so the code will not run outside of the computing environment provided by [Statistics Denmark](https://www.dst.dk/en/)

- [`1_format_data.R`](R/1_format_data.R): 
  This program creates datasets for analyzing the four aspects of disease progression. 
  Requirements: A `data.table` object named `data` holding all data for the project.
  
- [`2_fine_gray.R`](R/2_fine_gray.R): 
  This program runs analyses for diagnostic recurrence and pairwise diagnostic shifts. 
  Requirements: The `data.table` objects `data_recurrence` and `data_pairwise_shift` created in `1_format_data.R`
  
- [`3_ghosh_lin.R`](R/3_ghosh_lin.R): 
  This program runs analyses for the number of diagnostic shifts and the number of psychiatric contacts. 
  Requirements: The `data.table` objects `data_n_contacts` and `data_n_shifts` created in `1_format_data.R`




## Requirements

- [R](https://www.r-project.org/)
- [RStudio](https://posit.co/download/rstudio-desktop/) (optional but recommended)
- Required R packages all available on [CRAN](https://cran.r-project.org/) are explicitly loaded in each of the R-scripts


### Dataset

Most variables have a prefix indicating the type of variable:

* `d_`: date (class 'Date')
* `f_`: factor
* `i_`: integer

The programs require a dataset named `data` as a `data.table` object with one row for each combination of individual and contact and the following variables:

* `id`: Identification number of the individual 
* `d_in`: Date of first psychiatric contact 
* `d_exit`: Date of emigration, administrative censoring or death 
* `i_death`: Binary indicator variable denoting whether a death occurred at time d_exit 
* `f_first_diagnosis`: Diagnostic group at first contact
* `f_sex`: Sex of the individual
* `f_age`: Age group of the individual at first contact
* `f_calendar_time`: Calendar time period at first contact
* `d_start`: Start date of current psychiatric contact
* `f_diagnosis`: Diagnostic group representing the current contact

