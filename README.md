# A NEW BAYESIAN ECOLOGICAL INFERENCE MODEL WITH COVARIATES
Proposal of a new ecological inference method and comparison with previous ecological inference techniques.
We propose a new Bayesian ecological inference method to model the behaviour of individuals from the observed aggregated measures of the population with the help of covariates.
Additionally, we compare its results with the ecological inference methods from the current literature.

## How to use

### Requirements

* R, with the following packages:
  * glue
  * sf
  * tidyr
  * spdep 
  * lphom
  * ecolRxC
  * ei.Datasets
* R Studio
* JAGS
* WinBUGS
* Python
* Jupyter

### Execution

To generate the comparisons, execute in R Studio either `src/louisiana.r`, `src/new_zealand.r`.
To generate a latex version of the resulting tables, execute in R Studio `src/generate_table.r`.
To generate the plots from the resulting tables, execute in jupyter `src/Plotting.ipynb`.

## Folder structure
* Source code, `src/`:
  * Run code for a dataset:
    * `src/louisiana.r`: Execute the Louisiana comparisions for all models (~1 hour)
    * `src/new_zealand.r`: Execute the New Zealand comparisons for one of the two electorates (3~4.5 hours)
  * Code to run models:
    * `src/covariate.r`: Covariate Bayesian method
    * `src/spatial.r`: Spatial Bayesian method
    * `src/cluster.r`: Cluster Bayesian method
    * `src/wrapper.r`: Wrapper for linear programming method and latent structure method
    * `src/logit_covariate.r`: Our covariate Bayesian method, both with and without covariates
  * Auxiliary code:
    * `src/utilities.r`: Helper function to reduce code duplication.
    * `src/generate_table.r`: Latex tables for showcasing results
  * Plotting:
    * `src/plotting.r`: Simple plots for testing
    * `src/Plotting.ipynb`: Nice plots for showcasing results
  * Output, `src/out/`:
    * Louisiana results: `src/out/louisiana_stats.csv`, `src/out/louisiana.RData`
    * New Zealand, Auckland Central results: `src/out/New Zealand-Auckland Central_stats.csv`, `src/out/New Zealand-Auckland Central.RData`
    * New Zealand, Waiariki results: `src/out/New Zealand-Waiariki._stats.csv`, `src/out/New Zealand-Waiariki.RData`
* Dataset files, `data/`:
  * Louisiana dataset:
    * `data/US-counties.zip`: Map of the US, autmoatically extracted by the program
    * `data/Louisiana-2020_presidential_election.csv.zip`: Results of the 2020 preidential elections disagregated by party, race and parish
    * `data/Louisiana-2020_income_adjusted.csv`: Income statistics for Louisiana in 2020 disagregated by parish
    * `data/Louisiana-2020_education.csv`: Education statistics for Louisiana in 2020 disagregated by parish and demographic
  * New Zealand dataset is available in the `ei.Datasets` R package
* `models/`
  * Bayesian models: 
    * Covariate method: `models/covariate.txt`, `models/covariate-no_args.txt`
    * Spatial method: `models/covariate.txt`, `models/covariate-null_effect.txt`
    * Cluster method: `models/cluster.txt`, `models/cluster-1_cluster.txt`
    * Our covariate method: `models/logit_covariate.txt`, `models/logit_covariate-no_args.txt`, `models/logit_covariate-no_effects.txt`, `models/logit_covariate-no_effects-no_args.txt`
  * Simulation:
    * CAR-normal: `models/car_normal.txt`
    * Generate artificial data: `models/generate-local.txt`
