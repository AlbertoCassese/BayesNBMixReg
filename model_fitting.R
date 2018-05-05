# ***README***
# The following script is used to fit RNA-Seq and DNA methylation data to the Bayesian
# negative binomial mixture regression model proposed in the submitted manuscript titled 
# "Bayesian Negative Binomial Mixture Regression Models for the Analysis of Sequence 
# Count and Methylation Data."

# Before running the following code, please first load the data using "data_loader.R". The
# necessary inputs should be 
# (1) a n-by-p count matrix, i.e. RNA-Seq data, denoted by Y
# (2) a n-by-R matrix, i.e. DNA methylation data, denoted by X
# (3) a n-dimensional vector, i.e. sample allocation, denoted by z, note that z should 
# be a natural number starting at 0
# (4) a p-by-p matrix, i.e. gene-gene network, denoted by G
# (5) a R-by-p matrix, which restricts covariate selection, denoted by Theta, note that
# all of the elements in the matrix should be one if there is no such a restriction.
# ***END***

# Load library
require(Rcpp);

# Load function
Rcpp::sourceCpp('functions_x.cpp');



# ========================================================================================
# ========================================================================================
# Load data
# ========================================================================================
# ========================================================================================
# load("real_data_sub.RData");
# load("real_data_full.RData");
# load("simulated_data_test.RData"); # The true value of model parameters is also provided
load("synthetic_data_test.RData"); # The true value of model parameters is also provided



# ========================================================================================
# ========================================================================================
# Load algorithm setting
# ========================================================================================
# ========================================================================================
iter <- 10000;
burn <- iter/2;



# ========================================================================================
# ========================================================================================
# Load hyperparameters
# ========================================================================================
# ========================================================================================
# MRF prior for gamma
d <- -4;
f <- 1;

# Prior for omega
a_omega <- 0.2;
b_omega <- 1.8;

# Prior for phi
a_phi <- 1;
b_phi <- 0.01;

# Prior for sigma_0, sigma_alpha, sigma_beta
a_sigma <- 2.0;
b_sigma <- 1.0;



# ========================================================================================
# ========================================================================================
# Load initial configuration
# ========================================================================================
# ========================================================================================
n <- dim(Y)[1];
p <- dim(Y)[2];
R <- dim(X)[2];
K <- length(unique(z));
gamma_start <- rbinom(p, 1, 0.1);
phi_start <- rgamma(p, 10, 1);
alpha_0_start <- rnorm(p, 0, 1);
Alpha_start <- matrix(0, nrow = K, ncol = p);
Alpha_start[2:K, which(gamma_start == 1)] <- rnorm(sum(gamma_start == 1)*(K - 1), 0, 0.1);
Delta_start <- matrix(0L, nrow = R, ncol = p);
for (j in 1:p) {
  Delta_start[which(Theta[, j] == 1), j] <- rbinom(sum(Theta[, j] == 1), 1, 0.1);
}
Beta_start <- matrix(0, nrow = R, ncol = p);
Beta_start[which(Delta_start == 1)] <- rnorm(sum(Delta_start == 1), 0, 1);



# ========================================================================================
# ========================================================================================
# Implement MCMC algorithm
# ========================================================================================
# ========================================================================================
start_time <- proc.time();
M <- mcmc(iter, burn, d, f, a_omega, b_omega, a_phi, b_phi, a_sigma, b_sigma, X, Y, z, Theta, G, gamma_start, phi_start, alpha_0_start, Alpha_start, Delta_start, Beta_start);
end_time <- proc.time();
time <- end_time - start_time;
# The MCMC outputs are stored in M
# $phi:       the mean of MCMC samples of phi over iterations after burn in
# $alpha_0:   the mean of MCMC samples of alpha_0 over iterations after burn in
# $Alpha:     the mean of MCMC samples of Alpha over iterations after burn in
# $Beta:      the mean of MCMC samples of Beta over iterations after burn in
# $accept:    the acceptance rate for each updates: phi, alpha_0, (gamma, Alpha), Alpha, (Delta, Beta), Beta
# $gamma_sum: the total number of selected genes
# $Delta_sum: the total number of selected covariates for different genes
# $gamma_ppi: the marginal posterior probability of inclusion (PPI) of gamma
# $Delta_ppi: the marginal posterior probability of inclusion (PPI) of Delta



# ========================================================================================
# ========================================================================================
# Summarize result
# ========================================================================================
# ========================================================================================
print(paste0("Runtime = ", round(time[3], 1), "s"));
# print(paste0("Acceptance rate (phi): ", round(M$accept[1], 3)));
# print(paste0("Acceptance rate (alpha_0): ", round(M$accept[2], 3)));
# print(paste0("Acceptance rate (gamma, Alpha): ", round(M$accept[3], 3)));
# print(paste0("Acceptance rate (Alpha): ", round(M$accept[4], 3)));
# print(paste0("Acceptance rate (Delta, Beta): ", round(M$accept[5], 3)));
# print(paste0("Acceptance rate (Beta): ", round(M$accept[6], 3)));

plot(M$gamma_sum, type = "l", xlab = "Iterations", ylab = "Number of selected variables (genes)");
plot(M$gamma_ppi, type = "h", xlab = "Variable (gene) index", ylab = "Posterior probability of inclusion");

plot(M$Delta_sum, type = "l", xlab = "Iterations", ylab = "Number of selected covariates (probes)");
