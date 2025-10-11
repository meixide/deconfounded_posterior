# Predictive Posteriors under Hidden Confounding

A Bayesian framework for robust prediction across external domains in the presence of hidden confounders. This repository contains the implementation and reproducible code for our method, which provides well-calibrated predictive distributions with principled uncertainty quantification.
By Carlos García Meixide and David Ríos Insua

## Overview

Predicting outcomes in external domains remains a fundamental challenge when hidden confounders influence both predictors and outcomes. Traditional methods often require stringent assumptions, prior knowledge of distribution shifts, or rely on regularization strategies that compromise both estimation and predictive accuracy.

Our approach addresses these limitations by introducing a Bayesian framework that:

- **Generates well-calibrated predictive distributions** across unseen domains
- **Supports valid model inference** 
- **Achieves improved posterior contraction rates** as the number of observed datasets increases
- **Maintains empirical coverage** from low- to moderate-dimensional settings

## Repository Structure

### Core Implementation

- **`gi_pred.stan`** - One-dimensional, single-source version of BGI
- **`gi_hd.stan`** - Multivariate, multi-source version of BGI

### Analysis Scripts

- **`basic.R`** - Generates Figures 1, 2, and 3 from the paper. Ideal for exploring the algorithm and building intuition about the method.

- **`main_simulations.R`** - Reproduces the simulations from Section 5.2. Computationally intensive and optimized for parallel execution.

## Getting Started

### Prerequisites

- R (≥ 4.0)
- Stan and RStan
- Required R packages: (list your dependencies)

### Running the Examples

**Large-scale simulations:**

For HPC environments, launch parallel simulations with:
```bash
Rscript main_simulations.R <sample_size> <n_variables> <n_cores> <n_runs>
```

Example:
```bash
Rscript main_simulations.R 1000 10 8 24
```

Where:
- `1000` = sample size
- `10` = number of variables
- `8` = number of cores
- `24` = number of simulation runs

## Method Highlights


## Citation

If you use this code in your research, please cite:


## Contact

For questions or issues, please open an issue on GitHub or contact the authors.

---

*
