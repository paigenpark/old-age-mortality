# Old-Age Mortality Smoothing
## Project Overview

**Purpose:** This project evaluates various statisical methods for smoothing old-age mortality rates in the Human Mortality Database. 

**Models Tested:** Current methods explored include classic parametric models (Beard, Kannisto, Gompertz, Makeham, Log-quadratic), mixture models, and p-splines. Methods in development include Bayesian hierarchical models and machine learning apporaches. 

**Evaluation Metrics:** AIC is metric used for assessment currently. Planning to expand to cross-validation metrics like MSE/accuracy. 

## Repository Strucutre 

**Data:** Raw data are mostly downloaded input data files from mortality.org. Processed data include raw death data and exposure data calculated using the extinct cohort method. *More to come here, working on organizing data files.*

**Code:** Code files include data preparation, EDA, and function scripts. Also includes the current draft of the manuscript. 
*Data Preparation*
```
code/create_cohort_validation_data.qmd
```
code/prep_france.qmd
```
code/prep_scandi.qmd

*Exploratory Data Analysis*
```
code/eda.qmd

*Model fitting functions*
```
code/model_fitting_functions_cohort.R
```
code/local_flex_models.R

*Manuscript*
```
code/manuscript.qmd

**Figures:** Charts and tables generated during analysis.
