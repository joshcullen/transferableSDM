

## Simulate data (time interval) from generative model ##

# mu = dist * exp(b0_id + b1*bathym + b2*chla + b3*sst)
# a = mu * b

# dt ~ gamma(a, b)

set.seed(2022)

# Define data
ID <- rep(1:10, each = 1000)
N <- length(ID)

# Define coefficients
mu_bar <- rnorm(1, 0, 1)
sigma <- rexp(1, 1)
b0_id <- rnorm(n_distinct(ID), mean = mu_bar, sd = sigma)
b1 <- -0.5
b2 <- 0
b3 <- 0.8
b <- 0.2

# Generate time interval (dt) using subset of real data for `dist` and environ covars
for (i in 1:N) {
  mu <- dat3$dist * exp(b0_id[ID] + b1*dat3$bathym.s + b2*dat3$chla.s + b3*dat3$sst.s)
  a <- mu * b
  dt <- rgamma(N, shape = a, rate = b)
}



dat.list <- list(
  N = N,
  ID = ID,
  dt = dt,
  dist = dat3$dist,
  bathym = dat3$bathym,
  chla = dat3$chla,
  sst = dat3$sst
)

stan.model <- '
data {
  int N;                                  // sample size
  int ID[N];                           // ID label for each step
  vector[N] dt;                           // time interval (min)
  vector[N] dist;                         // Distance traveled for given step (m)
  vector[N] bathym;                       // Bathymetric depth (m)
  vector[N] chla;                         // Chlorophyll a concentration (mg m^-3 d^-1)
  vector[N] sst;                          // Sea surface temperature (C)
}


transformed data {
  vector[N] bathym_s;
  vector[N] chla_s;
  vector[N] sst_s;

  // Center and scale covariates
  bathym_s = (bathym - mean(bathym)) / sd(bathym);
  chla_s = (chla - mean(chla)) / sd(chla);
  sst_s = (sst - mean(sst)) / sd(sst);
}


parameters {
    real<lower=0> b;
    vector[10] b0_id;
    real b1;
    real b2;
    real b3;

    real mean1;
    real<lower=0> sd1;
}


transformed parameters {
  vector[N] mu;
  vector[N] a;

  for (i in 1:N){
    // mean of gamma distribution
    mu[i] = dist[i] * exp(b0_id[ID[i]] + b1*bathym_s[i] + b2*chla_s[i] + b3*sst_s[i]);
    // calculate the corresponding a[i] parameter
    a[i] = mu[i] * b;
  }
}


model {
  // likelihood
  for (i in 1:N){
    // b[i] = exp(g0+g1*Miss[i]);

    // likelihood
    dt[i] ~ gamma(a[i], b);
  }


  // priors
  for (j in 1:10){
    b0_id[j] ~ normal(mean1,sd1);
  }

  b1 ~ normal(0,1);
  b2 ~ normal(0,1);
  b3 ~ normal(0,1);

  b ~ exponential(1);

  mean1 ~ normal(0,1);
  sd1 ~ exponential(1);
}



'

# generated quantities {
#   vector[N] y_hat;
#
#   for (i in 1:N){
#     // mean of gamma distribution
#     y_hat[i] = dist[i] * exp(b0_id[ID[i]] + b1*bathym_s[i] + b2*chla_s[i] + b3*sst_s[i]);
#   }
# }

mod1 <- stan(model_code = stan.model, data = dat.list, chains = 4, iter = 2000, warmup = 1000, seed = 8675309)

params <- c('b','b0','b1','b2','b3')
print(mod1, digits_summary = 3, pars = params, probs = c(0.025, 0.5, 0.975))
