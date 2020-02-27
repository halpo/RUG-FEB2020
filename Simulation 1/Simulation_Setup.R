
# Description -----------------------------------------------------
# Simulate multivariate data then fit with variable selection
# Compute metrics for selection of variables and coverage of
# confidence intervals.
#



# Setup ------------------------------------------------------

library(elasticnet)
library(tidyverse)
library(assertthat)


# Functions -------------------------------------------------------

with_timing <- function(expr){
    time <- system.time(val <- force(expr))
    structure(val, time=time)
}
generate_data <-
function( n                 #< Number of observations
        , beta              #< Vector of Coefficients
        , Sigma             #< Variable correlation matrix
        ){
    if (is_formula(Sigma)) Sigma <- rlang::as_function(Sigma)() else
    if (is_function(Sigma)) Sigma <- Sigma
    assert_that( nrow(Sigma) == ncol(Sigma)
               , nrow(Sigma) == length(beta)
               , all(diag(Sigma) == 1)
               )


    # Generate raw data
	X <- mvtnorm::rmvnorm( n
	                     , mean = rep(0,length(beta))
	                     , sigma = Sigma
	                     )
    # Normalize for comparison
	X  <- scale(X,center=TRUE,scale=TRUE)

	# Compute response
	y  <- as.vector(X %*% beta)
	# add error
	ye <- y + rnorm(n,0,1)

	# Combind & Return
	structure( tibble( y
	                 , ye
                     , X
	                 )
	         , beta=beta
	         , true_response = y
	         )
}
analyze_with_enet<-function(data, beta, lambda=0){
	Y <-data$ye
	X  <-data$X
	# ps<-sum(beta!=0)
	# p<-length(beta)
	# n<-NROW(X)

	analysis <- enet(X, Y, lambda=lambda)
	cv     <- cv.enet( X, Y, K=5, lambda=lambda, s=seq(0,1,by=0.1)
	                 , plot.it=FALSE, mode="fraction", se=FALSE)
	minind <- which.min(cv$cv)
	smin   <- cv$s[minind]

	# fit  <- predict(analysis, X, s=smin,mode="fraction", type="fit")$fit
	coef <- predict(analysis, X, s=smin,mode="fraction", type="coef")$coef

	compute_performance_statistics(Y, X, coef, beta)
}
compute_performance_statistics<-function(Y,X,coef,beta){
	# Y the response vector
	# X the full set of predictors including the intercept in the first position.
	# coef the full set of estimated coefficients for the model including the in.
	# beta the full set of predictors including the intercept in the first position

	#validity checks
	assert_that( ncol(X) == length(coef)
	           , length(coef) == length(beta)
    	       , NROW(Y) == NROW(X)
    	       )
	# if(any(X[,1]!=1))stop("X must include the intercept")

	#false positives and false negatives
	false_positives=sum((coef!=0)&(beta==0))/sum(beta==0)
	false_negatives=sum((coef==0)&(beta!=0))/sum(beta!=0)

	X_s=X[,coef!=0,drop=FALSE]

	if(!all(coef==0)){
	# Estimates
	# yields the results of the hypothesis test for comparing against zero and the true value for the two methods.
	# there is overlap for the beta=0 coeficients.
	# only computed for the coeficients where the estimated coefficient is non zero.

	#Algorithm 1
		beta_adapt		= coef[coef!=0]
		rss_adapt		= crossprod(Y-(X_s%*%beta_adapt))
		mse_adapt  		= as.numeric(rss_adapt/(length(Y)-length(beta_adapt)))
		sigma_adapt	 	= diag(mse_adapt*solve(crossprod(X_s)))
		sd_adapt		= sqrt(sigma_adapt)
		adapt_zerotests = pt(-abs(beta_adapt/sd_adapt),df=length(Y)-length(beta_adapt))<.025
		# if(coef[1]!=0)adapt_zerotests<-adapt_zerotests[-1]
		adapt_truetests = pt(-abs((beta_adapt-beta[coef!=0])/sd_adapt),df=length(Y)-length(beta_adapt))<.025
		# if(coef[1]!=0)adapt_truetests<-adapt_truetests[-1]
	#Algorithm 2
		beta_oracle	 = coef(lm(Y~X_s-1))
		rss_oracle	 = crossprod(Y-X_s%*%beta_oracle)
		mse_oracle	 = as.numeric(rss_oracle/(length(Y)-length(beta_oracle)))
		sigma_oracle = diag(mse_oracle*solve(crossprod(X_s)))
		sd_oracle	 = sqrt(sigma_oracle)
		oracle_zerotests = pt(-abs(beta_oracle/sd_oracle),df=length(Y)-length(beta_oracle))<.025
		# if(coef[1]!=0)oracle_zerotests<-oracle_zerotests[-1]
		oracle_truetests = pt(-abs((beta_oracle-beta[coef!=0])/sd_oracle),df=length(Y)-length(beta_oracle))<.025
		# if(coef[1]!=0)oracle_truetests<-oracle_truetests[-1]

	#filter out intercept.
	#from here on the coef, beta, and test vectors are for the predictors only
		coef<-coef[-1]
		beta<-beta[-1]

	#zero beta CI coverage of zero
		adapt_zero_beta_CI_coverage_of_zero = (sum(!adapt_zerotests[beta[coef!=0]==0])+sum(coef[beta==0]==0))/sum(beta==0)
		oracle_zero_beta_CI_coverage_of_zero = (sum(!oracle_zerotests[beta[coef!=0]==0])+sum(coef[beta==0]==0))/sum(beta==0)

	#non-zero beta CI coverage of zero
		adapt_non_zero_beta_CI_coverage_of_zero  = (sum(!adapt_zerotests[beta[coef!=0]!=0])+sum(coef[beta!=0]==0))/sum(beta!=0)
		oracle_non_zero_beta_CI_coverage_of_zero = (sum(!oracle_zerotests[beta[coef!=0]!=0])+sum(coef[beta!=0]==0))/sum(beta!=0)

	#non-zero beta CI coverage of true beta
		adapt_non_zero_CI_coverage_of_beta = sum(!adapt_truetests[beta[coef!=0]!=0])/sum(beta!=0)
		oracle_non_zero_CI_coverage_of_beta = sum(!oracle_truetests[beta[coef!=0]!=0])/sum(beta!=0)
	} else { #special case for when no variables are selected
		false_positives = 0
		false_negatives = sum((beta!=0))
		adapt_zero_beta_CI_coverage_of_zero      = 1
		oracle_zero_beta_CI_coverage_of_zero     = 1
		adapt_non_zero_beta_CI_coverage_of_zero  = 1
		oracle_non_zero_beta_CI_coverage_of_zero = 1
		adapt_non_zero_CI_coverage_of_beta       = 0
		oracle_non_zero_CI_coverage_of_beta      = 0
	}

	#return value
    tibble( "rate of F+" = false_positives
		  , "rate of F-" = false_negatives
		  , "&beta;=0 CI coverage of zero (adapt)" = adapt_zero_beta_CI_coverage_of_zero
		  , "&beta;=0 CI coverage of zero (oracle)" = oracle_zero_beta_CI_coverage_of_zero
		  , "&beta;&ne;0 CI coverage of zero (adapt)" = adapt_non_zero_beta_CI_coverage_of_zero
		  , "&beta;&ne;0 CI coverage of zero (oracle)" = oracle_non_zero_beta_CI_coverage_of_zero
		  , "&beta;&ne;0 CI coverage of true beta (adapt)" = adapt_non_zero_CI_coverage_of_beta
		  , "&beta;&ne;0 CI coverage of true beta (oracle)" = oracle_non_zero_CI_coverage_of_beta
		  )
}
simulate_one <- function(n, beta, Sigma, analyze, ...){
    if (is_formula(beta)) beta <- rlang::as_function(beta)
    if (is_function(beta)) beta <- beta()
    assert_that(is.numeric(beta))
    analyze <- rlang::as_function(analyze)

    data <- with_timing(generate_data(n, beta, Sigma))
    # analysis
    analysis <- with_timing(suppressWarnings(analyze(data, beta, ...)))

    tibble( n=n
          , data.generation.time = list(attr(data, 'time'))
          , analysis.time = list(attr(analysis, 'time'))
          , analysis=list(analysis)
          , ...)
}
run_simulation <-function(m,...){
	rerun(m, simulate_one(...))
}

# Run Simulation --------------------------------------------------

if(FALSE){ # First Simulation
    p  <- 50
    ps <- p/5
    magnitude <- 1

    beta  <- sample(magnitude*c(rep(c(-1,1),each=ps/2),c(rep(0,p-ps))))
    Sigma <- rho**abs(outer(seq_along(beta), seq_along(beta), `-`))

    data <- generate_data(250, beta, Sigma)

    system.time({
        results <- simulate_one( 250, beta, Sigma
                               , analyze=analyze_with_enet
                               )
    })

    system.time({
    results <- run_simulation( m=100, n=250
                             , beta=beta, Sigma=Sigma
                             , analyze=analyze_with_enet)
    })



    m=1000  # number of simulation runs to perform

    simgrid  <-
        expand.grid(n=c(250, 500),p=c(p1,p2))

    save(simresults,file=paste("Sim Results ", format(Sys.time(), "%d-%m-%Y(%H,%M,%S)"),".rda",sep=''))
}
