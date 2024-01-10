# On the testing of multiple hypothesis using multiple data splitting


## Installation
```R
install_github("BioAlgs/MMM")
```

## Example
We use the following example to illustrate the basic usage of the MMM package. The MMM method is designed for any design matrix. We give several examples based on design matrix with different correlation structures. For simplicity, we will use synthetic data constructed from a linear model such that the response only depends on a small fraction of the variables.

In this example, we generate the design matrix with i.i.d. rows, and eachrow is generated from AR(1) model with autocorrelation 0.2. 


```R
n <- 300  # Number of observations
p <- 1000  # Number of predictors included in model
q <- 60  # Number of true predictors
amplitude = 5 # signal amplitude (for noise level = 1) 


# Generate the variables from a multivariate normal distribution
a1_nonnull = rep(amplitude/sqrt(n),q)
a2_nonnull = rep(amplitude/sqrt(n),q)
true_index= sample(p,q)
true_index=1:q
a1  = numeric(p)
a2  = numeric(p)
a1[true_index]=a1_nonnull
a2[true_index]=a2_nonnull
amat = cbind(a1,a2)
times <- 1:p
sigma <- 1
rho = 0.2
H <- abs(outer(times, times, "-"))
V <- sigma * rho^H
x = mvrnorm(n, mu = rep(0,p), Sigma=V)

# Generate the response.
y = sin(x %*% amat[,1]) + (x%*% amat[,2])^3 + rnorm(n, sd = 1);

# Fit the Gaussian mirror model
fit = gm(y, x, ncores=4)
# fit = gm(y, x, ncores=4,var.est=T) # Run GM with estiamtion of the variance of FDR

```
