
### single data-splitting method
analys = function(mm, ww, q){
  t_set = max(ww)
  for(t in ww){
    ps = length(mm[mm>=t])
    ng = length(na.omit(mm[mm<=-t]))
    rto = (ng)/max(ps, 1)
    if(rto<=q){
      t_set = c(t_set, t)
    }
  }
  thre = min(t_set)
  nz_est = which(mm>thre)
  nz_est
}

DS <- function(X, y, q){

  beta1_mat = matrix(0, 100, p)
  beta2_mat = matrix(0, 100, p)

  for(k in 1:10){
  n = dim(X)[1]; p = dim(X)[2]
  while(T){
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    ### get penalty lambda
    cvfit <- cv.glmnet(X[sample_index1, ], y[sample_index1], type.measure = "mse", nfolds = 10)
    lambda <- cvfit$lambda.1se
    beta1 <- as.vector(glmnet(X[sample_index1, ], y[sample_index1], family = "gaussian", alpha = 1, lambda = 0)$beta)
    nonzero_index <- which(beta1 != 0)
    if(length(nonzero_index)!=0)
      break
  }
  beta2 <- rep(0, p)
  beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)
  
  ### calculate the test statistics
  M <- abs(beta1 + beta2) - abs(beta1 - beta2)
  beta1_mat[k,] = beta1
  beta2_mat[k,] = beta2
  }

  M <- abs(beta1 + beta2) - abs(beta1 - beta2)
  beta1 = apply(beta1_mat, 2, mean)
  beta2 = apply(beta2_mat, 2, mean)
  
  selected_index <- analys(M, abs(M), q)
  
  # M_avg = apply(M_mat, 2, mean)
  return(list(selected_index = selected_index, M_avg= M))
}

library(MASS)
library(glmnet)
set.seed(123456)  # Set seed for reproducibility
n <- 1000  # Number of observations
p <- 300  # Number of predictors included in model
q <- 60  # Number of true predictors


nrep = 1
w_mat = matrix(0, p, nrep)
# rho_vec = seq(0,.4,.2)
rho_vec=0
amplitude=20

pow_ds = fdr_ds = matrix(0,nrep, length(rho_vec))

for(a in 1:length(rho_vec)){
for(r in 1: nrep){
    times <- 1:p
    sigma <- 1
    rho = rho_vec[a]
    V= matrix(rho, p,p)
    diag(V)= 1
    x = mvrnorm(n, mu = rep(0,p), Sigma=V)
    beta_q = rnorm(q, 0, amplitude)/sqrt(n)
    true_index= sample(p,q)
    beta = numeric(p)
    beta[true_index]=beta_q
    mu = x %*%  beta
    y = mu + rnorm(n)

    res=DS(x,y,0.1)

    selected.ds = unlist(res$selected_index)

    pow_ds[r,a] = length(intersect(selected.ds,true_index))/length(true_index)
	  fdr_ds[r,a] = 1-(length(intersect(selected.ds,true_index)))/length(selected.ds)
}
}

hist(res$M_avg[-true_index])