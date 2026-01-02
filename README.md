# Robust Genomic Prediction and Heritability Estimation using Density Power Divergence

This repository contains the complete implementation supporting the paper:

**Robust Genomic Prediction and Heritability Estimation using Density Power Divergence**  
Crop Science (2024)  
https://doi.org/10.1002/csc2.21430

The code provides robust and non-robust genomic prediction and heritability estimation methods based on one-stage and two-stage linear mixed models, with robustness achieved through Density Power Divergence (DPD).

---

## 📌 Overview

This repository includes:
- Simulation of genomic data under clean and contaminated scenarios
- One-stage and two-stage **non-robust** genomic prediction methods
- One-stage and two-stage **robust MDPDE-based** methods
- Application to real maize datasets
- Evaluation of prediction accuracy and heritability under outliers and model misspecification

---

## 📂 Repository Structure

data/
├── real/ # Real maize datasets
└── simulated/ # Generated via simulation scripts

simulation/
└── simulation_with_contamination.R

methods/
├── non_robust/
│ ├── simulated_one_two_stage.R
│ └── real_one_two_stage.R
└── mdpde/
├── simulated_mdpde_one_two_stage.R
└── real_mdpde_one_two_stage.R

results/ # Output placeholders (figures, tables)
environment/ # Session and package information


---

## ▶️ How to Run the Code

### 1. Simulate data (with and without contamination)
```r
source("simulation/simulation_with_contamination.R")
### Run non-robust methods
source("methods/non_robust/simulated_one_two_stage.R")
source("methods/non_robust/real_one_two_stage.R")
### Run proposed MDPDE-based robust methods
source("methods/mdpde/simulated_mdpde_one_two_stage.R")
source("methods/mdpde/real_mdpde_one_two_stage.R")
