# LongCOVID_Network

## Overview
This project performs **symptom network analysis of Long COVID** using co-occurrence-based network methods. It compares symptom networks between infected and uninfected groups, adjusting for confounders via inverse probability weighting (IPW).

## Analysis Pipeline

| Section | Description |
|---------|-------------|
| 0. Setup | Output directories and package loading |
| 1. Configuration | User-defined parameters (data path, symptom columns, covariates) |
| 2. Data preprocessing | Loading, derived group variables (age, BMI, CCI, vaccination) |
| 3. Helper functions | Binary coercion, symptom name cleaning, network utilities |
| 4A. Baseline network | Unweighted symptom co-occurrence network |
| 4B. IPW & Table 1 | Propensity score estimation and weighted descriptive statistics |
| 4C. IPW-weighted network | Confounding-adjusted symptom network |
| 4D. Time-bin analysis | Temporal symptom network across time periods |
| 4E. Bootstrap & permutation | Statistical inference (500 bootstrap / 1,000 permutation iterations) |

## Symptoms Analyzed
29 symptoms including fatigue, dyspnea, cognitive impairment, chest pain, palpitations, hair loss, sleep disturbance, loss of taste/smell, and others.

## Covariates for IPW
- Age group, Sex, BMI group
- Charlson Comorbidity Index (CCI)
- Smoking, Alcohol use
- Vaccination status, Marital status
- Household size, Education, Employment

## Requirements

### R Packages
```r
dplyr, tidyr, purrr, stringr, readxl,
igraph, ggplot2, ggrepel, ineq, tableone, survey
```

## Output
The pipeline generates:
- **Tables** — descriptive statistics, network centrality measures
- **Figures** — network visualizations, plots
- **Data** — processed datasets and network objects

## Usage
1. Update `data_path` in Section 1 to point to your data file
2. Verify symptom columns and covariate names match your dataset
3. Run the full script in R or RStudio
