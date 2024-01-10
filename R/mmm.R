#' On the testing of multiple hypothesis in sliced inverse regression
#'
#' 
#' @param y response variable. 
#' @param x input design matrix of dimension n x p; each row is an observation vector.
#' @param q target FDR level. The default value is 0.1.
#' @param nslices number of cores used for parallel computation of the mmm_1 to mmm_p. We recommend to use multiple cores which will greatly reduce the computational time.
#' @return mmm_selected: the index of selected variable.
#' @return mmm_statistics: the value for mirror statistics.

#' @author Xin Xing, \email{xin_xing@fas.harvard.edu} Zhigen Zhao, \email{zhaozhg@temple.edu}
#' @references \url{}
#' @keywords slice inverse regression, variable selection, fdr
#'
#' @examples
#' \donttest{n <- 300  # Number of observations
#'p <- 1000  # Number of predictors included in model
#'q <- 60  # Number of true predictors
#'amplitude = 5 # signal amplitude (for noise level = 1) 
#'
#'
#' # Generate the variables from a multivariate normal distribution
#' times <- 1:p
#'simmma <- 1
#' rho = 0.2
#'H <- abs(outer(times, times, "-"))
#'V <- simmma * rho^H
#'x = mvrnorm(n, mu = rep(0,p), Simmma=V)
#'
#'# Generate the response.
#'
#'a1_nonnull = rep(amplitude/sqrt(n),q)
#'a2_nonnull = rep(amplitude/sqrt(n),q)
#'true_index= sample(p,q)
#'true_index=1:q
#'a1  = numeric(p)
#'a2  = numeric(p)
#'a1[true_index]=a1_nonnull
#'a2[true_index]=a2_nonnull
#' y = sin(x %*% amat[,1]) + (x%*% amat[,2])^3 + rnorm(n, sd = 1);
#'
#'# Fit the MMM model
#'fit = mmm(y, x)
#'}
#'
#'
#' @export mmm
#' @importFrom MASS mvrnorm
#' @import dr



mmm = function(y, x, q=0.1, nslices= 8){

    n = length(y)
    p = dim(x)[2]
    x = scale(x)

	ind1 = sample(1:n, floor(n/2), replace=F)
	ind2 = setdiff(1:n, ind1)
	x1 = x[ind1, ]
	y1 = y[ind1]
	x2 = x[ind2,]
	y2 = y[ind2]

	fit1 = summary(dr(y1~x1[,],slice.function=dr.slices.arc,nslices=nslices, chi2approx="wood",numdir=2,method="sir"))
	fit2 = summary(dr(y2~x2[,],slice.function=dr.slices.arc,nslices=nslices, chi2approx="wood",numdir=2,method="sir"))


	proj1 = fit1$evectors %*% solve(t(fit1$evectors) %*% fit1$evectors ) %*% t(fit1$evectors)
	proj2 = fit2$evectors %*% solve(t(fit2$evectors) %*% fit2$evectors ) %*% t(fit2$evectors)


	mirror = numeric(p)
	for(i in 1:p){
		e = numeric(p)
		e[i] = 1
		e_proj1 = proj1 %*% e
		e_proj2 = proj2 %*% e
		# cos12 = sum(e_proj1 * e_proj2)/ (sqrt(sum(e_proj1^2)) * sqrt(sum(e_proj2^2)))
		cos12 = sum(e_proj1 * e_proj2)
		mirror[i] = cos12
	}


	w_vec_nonzero = mirror[which(mirror!=0)]
	w_idx = order(w_vec_nonzero, decreasing =T)
	est_fdp = numeric(length(w_vec_nonzero))
	true_fdp = numeric(length(w_vec_nonzero))
	for( i in 1:length(which(mirror>0))){
	    cut = w_vec_nonzero[w_idx[i]]
	    est_fdp[i] = sum(w_vec_nonzero < -cut)/(i+1)
	    true_fdp[i] = (i  - length(intersect(which(mirror!=0)[w_idx[1:i]], true_index)))/i
	}



	if(length(which(est_fdp-0.1>0))==0){
	    mmm_idx= length(which(mirror>0))
	}else{
	    mmm_idx = which(est_fdp-0.1>0)[1]-1
	}


	mmm_selected = which(mirror!=0)[w_idx[1:mmm_idx]]

return(results= list(mmm_selected=mmm_selected, mmm_statistics= mirror))



}

