# **trandiagnostic-patterns-in-mental disorders**

---

Paper: Nationwide transdiagnostic age and sex patterns in mental disorder disease developments

* Authors: [Eva N. S. Wandall](https://github.com/evaninawandall), [Per K. Andersen](https://researchprofiles.ku.dk/da/persons/per-kragh-andersen/), [Rune H. B. Christensen](https://github.com/runehaubo), [Michael E. Benros](https://researchprofiles.ku.dk/da/persons/jrd819-jrd819/)




### **R-scripts**

---

The R-scripts listed below are available in the R folder in this repository. 
The files contain code that reproduces key analyses of the paper. The data to 
run the code cannot be made publicly available so the code will not run outside 
of the computing environment provided by [Statistics Denmark](https://www.dst.dk/en/).

* `1_format_data.R`: This program creates data set for analysing the four aspects of disease progression. Requirements: data
* `2_fine_gray.R`: This program runs analyses for diagnostic recurrence and pairwise diagnostic shifts. Requirements: data_recurrence and data_pairwise_shift created in `1_format_data.R`
* `3_ghosh_lin.R`: This program runs analyses for the number of diagnostic shifts and the number of psychiatric contacts. Requirements: data_n_contacts and data_n_shifts created in `1_format_data.R`




### **Requirements**

---

* [R](https://www.r-project.org/)
* [RStudio](https://posit.co/download/rstudio-desktop/) (optional but recommended)
* Required R packages all available on [CRAN](https://cran.r-project.org/) are explicitly loaded in each of the R-scripts


##### **Dataset:**

Most variables have a prefix indicating the type of variable:

* `d_`: date (class 'Date')
* `f_`: factor
* `i_`: integer

The programs require the following dataset:

data: individual level data with the following variables (one record per individual per contact)

* `id`: Identification number of the individual
* `d_in`: Date of first psychiatric contact
* `d_exit`: Date of emigration, administrative censoring or death
* `i_death`: Binary indicator variable denoting whether a death occured at time d_exit
* `f_first_diagnosis`: Diagnostic group at first contact
* `f_sex`: Sex of the individual
* `f_age`: Age group of the individual at first contact
* `f_calendar_time`: Calendar time period at first contact
* `d_start`: Start date of current psychiatric contact
* `f_diagnosis`: Diagnostic group representing the current contact

