#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
NumericVector csample_num(NumericVector x, int size, bool replace);
NumericVector getIndex1(NumericVector beta);
NumericVector getIndex0(NumericVector beta, IntegerVector theta);
NumericVector getIndex(IntegerVector theta);
int getNum1(NumericVector beta);
int getNum0(NumericVector beta, IntegerVector theta);
NumericVector getValue1(NumericVector beta);
arma::mat Sigma(NumericMatrix X, NumericVector index, double h);

// [[Rcpp::export]]
Rcpp::List mcmc(int iter, int burn, double d, double f, double a_omega, double b_omega, double a_phi, double b_phi, double a_sigma, double b_sigma, NumericMatrix X, NumericMatrix Y, IntegerVector z, IntegerMatrix Theta, IntegerMatrix G, IntegerVector gamma_start, NumericVector phi_start, NumericVector alpha_0_start, NumericMatrix Alpha_start, IntegerMatrix Delta_start, NumericMatrix Beta_start) {
  int n = Y.nrow();
  int p = Y.ncol();
  int R = X.ncol();
  int K = max(z) + 1;
  double tau2 = 1;
  double ee = 0.1;
  int count = 1;
  int b, i, j, jj, k, r, e;
  double hastings, lambda, lambda_temp, phi_temp, alpha_0_temp, alpha_temp, beta_temp;
  IntegerVector gamma_1(p, 1);
  NumericVector alpha_array_temp(K);
  NumericVector beta_array_temp;
  NumericVector delta_index;
  NumericVector E(p);
  for (j = 0; j < p; j++)
  {
    for (r = 0; r < R; r++)
    {
      E(j) = E(j) + Theta(r, j);
    }
    E(j) = E(j)*0.1;
    if (E(j) <= 1)
    {
      E(j) = 1;
    }
    else
    {
      E(j) = floor(E(j));
    }
  }
  
  NumericVector propose(6);
  NumericVector accept(6);
  IntegerVector gamma_sum(iter, 0);
  IntegerVector Delta_sum(iter, 0);
  NumericVector gamma_ppi(p);
  NumericMatrix Delta_ppi(R, p);
  NumericVector phi_mean(p);
  NumericVector alpha_0_mean(p);
  NumericMatrix Alpha_mean(K, p);
  NumericMatrix Beta_mean(R, p);
  for (j = 0; j < p; j++)
  {
    phi_mean(j) = 0;
    alpha_0_mean(j) = 0;
    gamma_ppi(j) = 0;
    for (k = 0; k < K; k++)
    {
      Alpha_mean(k, j) = 0;
    }
    for (r = 0; r < R; r++)
    {
      Beta_mean(r, j) = 0;
      Delta_ppi(r, j) = 0;
    }
  }
  for (i = 0; i < 6; i++)
  {
    accept(i) = 0;
    propose(i) = 0;
  }
  
  // ================================================================================
  // MCMC Initialization ============================================================
  // ================================================================================
  NumericVector s(n, 1.0);
  IntegerVector gamma(p); 
  NumericVector phi(p);
  NumericVector alpha_0(p);
  NumericMatrix Alpha(K, p); 
  IntegerMatrix Delta(R, p); 
  NumericMatrix Beta(R, p);
  NumericMatrix reg(n, p);
  for (j = 0; j < p; j++)
  {
    gamma(j) = gamma_start(j);
    phi(j) = phi_start(j);
    alpha_0(j) = alpha_0_start(j);
    for (k = 0; k < K; k++)
    {
      Alpha(k, j) = Alpha_start(k, j);
    }
    for (r = 0; r < R; r++)
    {
      Delta(r, j) = Delta_start(r, j);
      Beta(r, j) = Beta_start(r, j);
    }
  }
  for (j = 0; j < p; j++)
  {
    for (i = 0; i < n; i++)
    {
      reg(i, j) = 0;
      for (r = 0; r < R; r++)
      {
        if(Delta(r, j) == 1) 
        {
          reg(i, j) = reg(i, j) + X(i, r)*Beta(r, j);
        }
      }
    }
  }
  
  // ================================================================================
  // MCMC algorithm =================================================================
  // ================================================================================
  for (b = 0; b < iter; b++) 
  {
    // Update phi -------------------------------------------------------------------
    for (j = 0; j < p; j++)
    {
      phi_temp = rgamma(1, a_phi, 1/b_phi)(0);
      hastings = 0;
      for (i = 0; i < n; i++)
      {
        lambda = s(i)*exp(alpha_0(j) + Alpha(z(i), j) + reg(i, j));
        hastings = hastings + (lgamma(Y(i, j) + phi_temp) - lgamma(phi_temp) + phi_temp*log(phi_temp) - (phi_temp + Y(i, j))*log(phi_temp + lambda));
        hastings = hastings - (lgamma(Y(i, j) + phi(j)) - lgamma(phi(j)) + phi(j)*log(phi(j)) - (phi(j) + Y(i, j))*log(phi(j) + lambda));
      }
      propose(0)++;
      if (hastings >= log(double(rand()%10001)/10000))
      {
        phi(j) = phi_temp;
        accept(0)++;
      }
    }
    
    // Update alpha_0 --------------------------------------------------------------
    for (j = 0; j < p; j++)
    {
      alpha_0_temp = rnorm(1, alpha_0(j), sqrt(ee*tau2))(0);
      hastings = 0;
      for (i = 0; i < n; i++)
      {
        lambda_temp = s(i)*exp(alpha_0_temp + Alpha(z(i), j) + reg(i, j));
        lambda = s(i)*exp(alpha_0(j) + Alpha(z(i), j) + reg(i, j));
        hastings = hastings + (Y(i, j)*alpha_0_temp - (phi(j) + Y(i, j))*log(phi(j) + lambda_temp));
        hastings = hastings - (Y(i, j)*alpha_0(j) - (phi(j) + Y(i, j))*log(phi(j) + lambda));
      }
      hastings = hastings + (-a_sigma - 1/2)*log(b_sigma + alpha_0_temp*alpha_0_temp/2);
      hastings = hastings - (-a_sigma - 1/2)*log(b_sigma + alpha_0(j)*alpha_0(j)/2);
      propose(1)++;
      if (hastings >= log(double(rand()%10001)/10000))
      {
        alpha_0(j) = alpha_0_temp;
        accept(1)++;
      }
    }
    
    // Update gamma and alpha ------------------------------------------------------
    for (e = 0; e < 5; e++)
    {
      j = rand()%p;
      hastings = 0;
      if (Alpha(1, j) == 0) // Add
      {
        alpha_array_temp = rnorm(K, 0, sqrt(tau2));
        alpha_array_temp(0) = 0;
        hastings = hastings + d;
        for (jj = 0; jj < p; jj++)
        {
          hastings = hastings + f*G(j, jj);
        }
        for (k = 1; k < K; k++)
        {
          hastings = hastings - (-0.5*log(2*M_PI) - 0.5*log(tau2) - alpha_array_temp(k)*alpha_array_temp(k)/(2*tau2));
          hastings = hastings + (-0.5*log(2*M_PI) + a_sigma*log(b_sigma) - lgamma(a_sigma) + lgamma(a_sigma + 0.5) - (a_sigma + 0.5)*log(b_sigma + alpha_array_temp(k)*alpha_array_temp(k)/2));
        }
      }
      else // Delete
      {
        alpha_array_temp(0) = 0;
        hastings = hastings - d;
        for (jj = 0; jj < p; jj++)
        {
          hastings = hastings - f*G(j, jj);
        }
        for (k = 1; k < K; k++)
        {
          alpha_array_temp(k) = 0;
          hastings = hastings + (-0.5*log(2*M_PI) - 0.5*log(tau2) - Alpha(k, j)*Alpha(k, j)/(2*tau2));
          hastings = hastings - (-0.5*log(2*M_PI) + a_sigma*log(b_sigma) - lgamma(a_sigma) + lgamma(a_sigma + 0.5) - (a_sigma + 0.5)*log(b_sigma + Alpha(k, j)*Alpha(k, j)/2));
        }
      }
      for (i = 0; i < n; i++)
      {
        lambda_temp = s(i)*exp(alpha_0(j) + alpha_array_temp(z(i)) + reg(i, j));
        lambda = s(i)*exp(alpha_0(j) + Alpha(z(i), j) + reg(i, j));
        hastings = hastings + (Y(i, j)*alpha_array_temp(z(i)) - (phi(j) + Y(i, j))*log(phi(j) + lambda_temp));
        hastings = hastings - (Y(i, j)*Alpha(z(i), j) - (phi(j) + Y(i, j))*log(phi(j) + lambda));
      }
      propose(2)++;
      if (hastings >= log(double(rand()%10001)/10000))
      {
        for (k = 0; k < K; k++)
        {
          Alpha(k, j) = alpha_array_temp(k);
        }
        accept(2)++;
      }
    }
    
    // Update alpha ----------------------------------------------------------------
    for (j = 0; j < p; j++)
    {
      if (Alpha(1, j) != 0)
      {
        for (k = 1; k < K; k++)
        {
          alpha_temp = rnorm(1, Alpha(k, j), sqrt(ee*tau2))(0);
          hastings = 0;
          for (i = 0; i < n; i++)
          {
            if (z(i) == k)
            {
              lambda_temp = s(i)*exp(alpha_0(j) + alpha_temp + reg(i, j));
              lambda = s(i)*exp(alpha_0(j) + Alpha(z(i), j) + reg(i, j));
              hastings = hastings + (Y(i, j)*alpha_temp - (phi(j) + Y(i, j))*log(phi(j) + lambda_temp));
              hastings = hastings - (Y(i, j)*Alpha(z(i), j) - (phi(j) + Y(i, j))*log(phi(j) + lambda));
            }
          }
          hastings = hastings + (-a_sigma - 1/2)*log(b_sigma + alpha_temp*alpha_temp/2);
          hastings = hastings - (-a_sigma - 1/2)*log(b_sigma + Alpha(k, j)*Alpha(k, j)/2);
          propose(3)++;
          if (hastings >= log(double(rand()%10001)/10000))
          {
            Alpha(k, j) = alpha_temp;
            accept(3)++;
          }
        }
      }
    }
    
    // Update delta and beta -------------------------------------------------------
    for (j = 0; j < p; j++)
    {
      for (e = 0; e < E(j); e++)
      {
        r = csample_num(getIndex(Theta(_, j)), 1, 0)(0);
        hastings = 0;
        if (Beta(r, j) == 0) // Add
        {
          beta_temp = rnorm(1, 0, sqrt(tau2))(0);
          hastings = hastings + (lgamma(a_omega + 1) + lgamma(b_omega));
          hastings = hastings - (lgamma(a_omega) + lgamma(b_omega + 1));
          hastings = hastings - (-0.5*log(2*M_PI) - 0.5*log(tau2) - beta_temp*beta_temp/(2*tau2));
          hastings = hastings + (-0.5*log(2*M_PI) + a_sigma*log(b_sigma) - lgamma(a_sigma) + lgamma(a_sigma + 0.5) - (a_sigma + 0.5)*log(b_sigma + beta_temp*beta_temp/2));
        }
        else // Delete
        {
          beta_temp = 0;
          hastings = hastings - (lgamma(a_omega + 1) + lgamma(b_omega));
          hastings = hastings + (lgamma(a_omega) + lgamma(b_omega + 1));
          hastings = hastings + (-0.5*log(2*M_PI) - 0.5*log(tau2) - Beta(r, j)*Beta(r, j)/(2*tau2));
          hastings = hastings - (-0.5*log(2*M_PI) + a_sigma*log(b_sigma) - lgamma(a_sigma) + lgamma(a_sigma + 0.5) - (a_sigma + 0.5)*log(b_sigma + Beta(r, j)*Beta(r, j)/2));
        }
        for (i = 0; i < n; i++)
        {
          lambda_temp = s(i)*exp(alpha_0(j) + Alpha(z(i), j) + reg(i, j) - X(i, r)*Beta(r, j) + X(i, r)*beta_temp);
          lambda = s(i)*exp(alpha_0(j) + Alpha(z(i), j) + reg(i, j));
          hastings = hastings + (Y(i, j)*(reg(i, j) - X(i, r)*Beta(r, j) + X(i, r)*beta_temp) - (phi(j) + Y(i, j))*log(phi(j) + lambda_temp));
          hastings = hastings - (Y(i, j)*reg(i, j) - (phi(j) + Y(i, j))*log(phi(j) + lambda));
        }
        propose(4)++;
        if (hastings >= log(double(rand()%10001)/10000))
        {
          for (i = 0; i < n; i++)
          {
            reg(i, j) = reg(i, j) - X(i, r)*Beta(r, j) + X(i, r)*beta_temp;
          }
          Beta(r, j) = beta_temp;
          accept(4)++;
        }
      }
    }
    
    // Update beta -----------------------------------------------------------------
    for (j = 0; j < p; j++)
    {
      // Update beta individually
      for (r = 0; r < R; r++)
      {
        if (Beta(r, j) != 0)
        {
          beta_temp = rnorm(1, Beta(r, j), sqrt(ee*tau2))(0);
          hastings = 0;
          for (i = 0; i < n; i++)
          {
            lambda_temp = s(i)*exp(alpha_0(j) + Alpha(z(i), j) + reg(i, j) - X(i, r)*Beta(r, j) + X(i, r)*beta_temp);
            lambda = s(i)*exp(alpha_0(j) + Alpha(z(i), j) + reg(i, j));
            hastings = hastings + (Y(i, j)*(reg(i, j) - X(i, r)*Beta(r, j) + X(i, r)*beta_temp) - (phi(j) + Y(i, j))*log(phi(j) + lambda_temp));
            hastings = hastings - (Y(i, j)*reg(i, j) - (phi(j) + Y(i, j))*log(phi(j) + lambda));
          }
          hastings = hastings + (-a_sigma - 1/2)*log(b_sigma + beta_temp*beta_temp/2);
          hastings = hastings - (-a_sigma - 1/2)*log(b_sigma + Beta(r, j)*Beta(r, j)/2);
          propose(5)++;
          if (hastings >= log(double(rand()%10001)/10000))
          {
            for (i = 0; i < n; i++)
            {
              reg(i, j) = reg(i, j) - X(i, r)*Beta(r, j) + X(i, r)*beta_temp;
            }
            Beta(r, j) = beta_temp;
            accept(5)++;
          }
        }
      }
    }
    
    // Monitor the process -----------------------------------------------------------------
    if ((b*100)/iter == count)
    {
      Rcout << count << "% has been done\n";
      count++;
    }
    if (b > burn)
    {
      for (j = 0; j < p; j++)
      {
        phi_mean(j) = phi_mean(j) + phi(j);
        alpha_0_mean(j) = alpha_0_mean(j) + alpha_0(j);
        for (k = 1; k < K; k++)
        {
          if (gamma(j) == 1)
          {
            Alpha_mean(k, j) = Alpha_mean(k, j) + Alpha(k, j);
          }
        }
        for (r = 0; r < R; r++)
        {
          if (Delta(r, j) == 1)
          {
            Beta_mean(r, j) = Beta_mean(r, j) + Beta(r, j);
          }
        }
      }
    }
    for (j = 0; j < p; j++)
    {
      if (Alpha(1, j) != 0)
      {
        gamma_sum(b)++;
        if (b > burn)
        {
          gamma_ppi(j)++;
        }
      }
      for (r = 0; r < R; r++)
      {
        if (Beta(r, j) != 0)
        {
          Delta_sum(b)++;
          if (b > burn)
          {
            Delta_ppi(r, j)++;
          }
        }
      }
    }
  }
  
  for (j = 0; j < p; j++)
  {
    phi_mean(j) = phi_mean(j)/(iter - burn);
    alpha_0_mean(j) = alpha_0_mean(j)/(iter - burn);
    gamma_ppi(j) = gamma_ppi(j)/(iter - burn);
    for (k = 1; k < K; k++)
    {
      Alpha_mean(k, j) = Alpha_mean(k, j)/(iter - burn);
    }
    for (r = 0; r < R; r++)
    {
      Beta_mean(r, j) = Beta_mean(r, j)/(iter - burn);
      Delta_ppi(r, j) = Delta_ppi(r, j)/(iter - burn);
    }
  }
  
  for (i = 0; i < 6; i++)
  {
    accept(i) = accept(i)/propose(i);
  }
  return Rcpp::List::create(Rcpp::Named("phi") = phi_mean, Rcpp::Named("alpha_0") = alpha_0_mean, Rcpp::Named("Alpha") = Alpha_mean, Rcpp::Named("Beta") = Beta_mean, Rcpp::Named("accept") = accept, Rcpp::Named("gamma_sum") = gamma_sum, Rcpp::Named("Delta_sum") = Delta_sum, Rcpp::Named("gamma_ppi") = gamma_ppi, Rcpp::Named("Delta_ppi") = Delta_ppi);
}

// [[Rcpp::export]]
NumericVector getIndex1(NumericVector beta) {
  int p = beta.size();
  int p_1 = 0;
  int count = 0;
  int j;
  for (j = 0; j < p; j++)
  {
    if (beta(j) != 0)
    {
      p_1++;
    }
  }
  NumericVector index(p_1);
  for (j = 0; j < p; j++)
  {
    if (beta(j) != 0)
    {
      index(count) = j;
      count++;
    }
  }
  return index;
}

// [[Rcpp::export]]
NumericVector getIndex0(NumericVector beta, IntegerVector theta) {
  int p = beta.size();
  int p_0 = 0;
  int count = 0;
  int j;
  for (j = 0; j < p; j++)
  {
    if (beta(j) == 0 && theta(j) == 1)
    {
      p_0++;
    }
  }
  NumericVector index(p_0);
  for (j = 0; j < p; j++)
  {
    if (beta(j) == 0 && theta(j) == 1)
    {
      index(count) = j;
      count++;
    }
  }
  return index;
}

// [[Rcpp::export]]
NumericVector getIndex(IntegerVector theta) {
  int p = theta.size();
  int p_1 = 0;
  int count = 0;
  int j;
  for (j = 0; j < p; j++)
  {
    if (theta(j) == 1)
    {
      p_1++;
    }
  }
  NumericVector index(p_1);
  for (j = 0; j < p; j++)
  {
    if (theta(j) == 1)
    {
      index(count) = j;
      count++;
    }
  }
  return index;
}

// [[Rcpp::export]]
int getNum1(NumericVector beta) {
  int p = beta.size();
  int p_1 = 0;
  for (int j = 0; j < p; j++)
  {
    if (beta(j) != 0)
    {
      p_1++;
    }
  }
  return p_1;
}

// [[Rcpp::export]]
int getNum0(NumericVector beta, IntegerVector theta) {
  int p = beta.size();
  int p_0 = 0;
  for (int j = 0; j < p; j++)
  {
    if (beta(j) == 0 && theta(j) == 1)
    {
      p_0++;
    }
  }
  return p_0;
}

// [[Rcpp::export]]
NumericVector getValue1(NumericVector beta) {
  int p = beta.size();
  int p_1 = 0;
  int count = 0;
  int j;
  for (j = 0; j < p; j++)
  {
    if (beta(j) != 0)
    {
      p_1++;
    }
  }
  NumericVector value(p_1);
  for (j = 0; j < p; j++)
  {
    if (beta(j) != 0)
    {
      value(count) = beta(j);
      count++;
    }
  }
  return value;
}

// [[Rcpp::export]]
NumericVector csample_num(NumericVector x, int size, bool replace) {
  NumericVector ret = RcppArmadillo::sample(x, size, replace);
  return ret;
}

//random multivariate normal sample generator using RcppArmadillo
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat Sigma) {
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y*arma::chol(Sigma);
}

// [[ Rcpp :: export ()]]
arma::mat Sigma(NumericMatrix X, NumericVector index, double h) {
  int n = X.nrow();
  int R = index.size();
  int count = 0;
  int i, r;
  arma::mat XX(n, R);
  arma::mat HH(R, R);
  for (r = 0; r < R; r++)
  {
    for (i = 0; i < n; i++)
    {
      XX(i, r) = X(i, index(count));
    }
    for (i = 0; i < R; i++)
    {
      HH(i, r) = h; 
    }
    count++;
  }
  return(HH%inv(XX.t()*XX));
}