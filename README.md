# BayesNBMixReg
This repository is used for accessing the performance of the Bayesian negative binomial mixture model, which was proposed in the submitted manuscript titled "Bayesian negative binomial mixture regression models for the analysis of sequence count and methylation data." Before running the code, please install R package Rcpp.

When running "model_fitting.R" to fit the proposed model, please first load the data. You can load your own data or one of the three examples in the repository.

*Note 1, the notations in the code and data follow the notations in the manuscript.

*Note 2, the results obtained by running the code in this repository may not be exactly the same as the results reported in the manuscript, because we reported the results by pooling multiple MCMC chains in the manuscript.

*Note 3, for large scale dataset, it takes longer time. Please first try small dataset such as "real_data_sub.RData", or run the MCMC algorithm with small number of iterations (e.g. 1,000).
