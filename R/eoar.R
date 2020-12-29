#' @export
#'
#' @title eoar - Evidence of Absence Regression model estimation.
#'
#' @description  An EoAR model
#' consists of a log-linear model for lambda, the mean number of search
#' targets per "cell", where
#' "cell" is a measured experimental unit such as a turbine, season, year, etc.
#' Inputs include
#' information about the number of targets found (\code{Y}) and the
#' g-values (=probatility of
#' discovery) at all
#' searched sites.  The method is Bayesian and allows either uniform
#' or informed priors for lambda.
#' Estimation is performed in JAGS.
#'
#' @param lambda A model formula for the lambda parameters of EoAR.
#' This formulat has the form \code{Y ~ X1 + X2} + etc. (exactly like \code{\link{lm}}).
#' Here,  \code{Y} is a vector containing number of carcasses found, one element per "cell".
#' For example, if
#' multiple sites are searched in a single season, elements of Y should be the number
#' of targets found at each site.  If multiple sites are searched during multiple seasons,
#' \code{Y} should contain the number of targets found at each site X season
#' combination.  Length
#' of \code{Y} is number of sites times number of seasons, assuming no missing.  Covariate
#' vectors (or matricies) \code{X1} etc. must have the same length as \code{Y} and
#' can be anything that \code{formula} handles (e.g., factors, interactions,
#' \code{I(x^2)}, etc.)
#'
#' @param beta.params A data frame containing, at a minimum, two columns named \code{$alpha}
#' and \code{$beta}.
#' These are the alpha and beta parameters of a Beta distribution that determine the g-values
#' in each "cell".
#' Length of \code{$alpha} and \code{$beta} is either one, in which case g is assumed constant
#' over "cells",
#' or one element per "cell".
#'
#' @param data An optional data frame, list or environment
#' (or object coercible by \code{as.data.frame} to a data frame)
#' containing all variables in the model. If variables are not found in data,
#' the variables are taken from environment(formula), typically the
#' environment from which \code{eoar} is called.
#'
#' @param offset An optional offset term(s) for the linear part of the model.
#' This should be NULL or a numeric vector of length equal to the number
#' of "cells". One or more offset terms can be included in the formula
#' (e.g., Y~ offset(log(a))+X), and if more than one is specified
#' their sum is used. See \code{\link{model.offset}}.  Due to the log link
#' used here, offsets should generally be logged first. When the log of an offset
#' is included, the linear part of the model can be read as
#' log(lambda/offset) ~ X1 + x2 + ...etc.
#'
#' @param priors An optional data frame specifying vague priors or
#' containing information
#' about (informed) prior distributions of coefficients in the lambda model.
#'
#' If \code{priors} = NULL, vague normal priors are assumed for
#' all coefficients.  Coefficients are on a log scale.  The
#' vague normal priors chosen here have mean 0 and very large standard deviations.
#' Check the output object to ensure the choosen standard deviations
#' are sufficiently large enough to be considered vague.
#'
#' If \code{priors} is not NULL, it should be a data frame containing
#' at a minimum \code{$mean} and \code{$sd} specifying the
#' mean and standard deviation of a normal prior for coefficients in the
#' lambda model.  Remember that coefficients are on the log scale, do not
#' specify means on the response (i.e., count) scale. Provide log(count) means.
#' The \code{row.names} of \code{priors} are matched to
#' coefficients, and this provides a way to use informed priors for some
#' coefficients and vague priors for others.  For example, if \code{X1} and
#' \code{X2} are both in the model and \code{priors}
#' contains only one row with \code{row.name} == "X1",  \code{$mean} and \code{$sd}
#' values for that row are used to inform the coefficient of \code{X1},
#' but vague priors are used for all other coefficients. See examples.
#'
#' @param conf.level The confidence level for all posterior confidence intervals.
#' Part of the value returned by this routine is \code{$intervals} containing
#' posterior confidence intervals for all parameters.  \code{conf.level} specifies
#' the two-tailed confidence level for all these confidence intervals.
#'
#' @param nburns The number of burn-in steps to take during MCMC sampling.
#'
#' @param niters The number of sampling steps to take during MCMC sampling.
#'
#' @param nthins The amount of MCMC chain thinning to do.  Every (\code{nthins})-th
#' sampling iteration in each chain is saved, while the
#' rest are discarded. Each chain has \code{niters} iterations. In
#' the end, a total of \code{nchains*niters/nthins} samples from
#' the posterior are available to the user.
#'
#' @param nchains The number of MCMC sampling chains. Must specify 2 or more to
#' check convergence.
#'
#' @param nadapt The number of adapting iterations to perform before
#' burn-in.  During adaptin, JAGS is trying to optimize it's proposal distribution
#' and stepsize to increase convergence speed.
#' 
#' @param computeIC A logical (default = FALSE) indicating whether to compute
#' and return information criteria. Currently, WAIC and the deviance of the
#' model are returned if \code{computeIC==TRUE}.
#'
#' @param quiet Logical indicating whether to print output during estimation.
#' \code{quiet==FALSE} shows text based progress bars during estimation.
#'
#' @param seeds A vector of length \code{nchains} containing random number
#' seeds for the MCMC sampler.  If NULL, \code{nchains} random numbers are
#' generated from R's base random number generator which is controled outside
#' this routine using \code{set.seed}.  Note
#' that \code{set.seed} has no effect on the random number sequences used
#' in JAGS because JAGS is a separate package.
#' The seeds, whether chosen by this routine or specified, are
#' stored in the output object.
#' 
#' @param doEoA A logical.
#'
#' @details
#' Observed quantities in the model are Y[i] = number of targets
#' observed in cell i, and the covariates X[1i], X[2i], etc. Parameters in
#' the model are M[i] = number of total targets in cell i including those
#' observed and those missed, lambda[i] = the mean number of (or rate of) targets
#' in cell i per offset[i], and g[i] = probability of observing a target in cell i.
#'
#' If a log(offset) term is not included, lambda[i] is mean number of targets per
#' cell.  For example, if a cell is one season of monitoring at single sites, lambda
#' is targets per season per site.  If a cell is one season of monitoring at multiple sites
#' (i.e., a facility), lambda is targets per season per facility.  If a log(offset) term
#' is included, lambda[i] is the cell-specific mean number of targets per offset unit.  For example, if
#' a cell contains counts from an entire  wind power facility and the offset is
#' log(number of turbines), lambda
#' is number of targets per turbine by facility in the analysis.
#' If the offset is log(days in season), lambda is
#' number of targets per day for each cell.  If two offset terms are included like
#' \code{offset(log(turbines)) + offset(log(days))}, lambda is number of targets
#' per turbine per day for the cell.
#'
#' Parameter \code{Mtot} is the sum of all \code{M} parameters in the analysis.
#' This derived parameters may not mean anything in specific problems, but is
#' a fairly common derived parameter.
#'
#' The EoAR model implemented here is :
#' \enumerate{
#'   \item Target count Y[i] is assumed to be binomial(M[i],g[i]).
#'   \item Binomial index M[i] is assumed to be poisson(lambda[i]).
#'   \item Binomial probability g[i] is assumed to be Beta(alpha[i],beta[i]).
#'   \item Poisson rate lambda[i] is constrained to be a log-linear function
#' of covariaes. That is, log(lambda[i]) = A[0] + A[1]x[1i] + A[2]x[2i] + etc. If
#' an log(offset) term is included, the model is log(lambda[i]/offset[i]) = A[0] + A[1]x[1i] + A[2]x[2i] + etc.
#'   \item Beta hyper-parameters alpha[i] and beta[i] are constants.
#' }
#'
#' There are two ways
#' to obtain the exact same results across multiple
#' runs.  One method is to specify
#' \code{seeds} here. This will set
#' the MCMC seeds in JAGS so that exact chains are reproduced.
#' For example, if \code{run1} is the
#' result of a previous call, \code{eoar(...,seeds=run1$seeds)}
#' will replicate \code{run1} exactly.
#' The second method is to use R's default \code{set.seed}
#' just before calling this routine.
#'
#' @return An object of class "eoar".  EoAR objects are lists containing the following components:
#' \itemize{
#'   \item \code{estimates} : a matrix containing parameter estimates and
#'   standard errors.  This matrix contains  one row per parameter and
#'   two columns.  \code{rownames(x$estimates)} lists the parameters.
#'   Columns are \code{Estimate}, containing the median value of
#'   of that parameter's posterior marginal, and \code{SD}, containing
#'   the standard deviation of the parameter's posterior marginal.  The
#'   following parameters are included:
#'   \itemize{
#'     \item \code{M} : Mortality estimates in individual "cells". For example,
#'     \code{M[4]} is the estimate (median of posterior) of mortalities at
#'     the site represented in cell 4 of the data set.
#'
#'     \item \code{Mtot} : A derived parameter, the median of the posterior
#'     distribution of the sum of all \code{M}'s.
#'
#'     \item \code{lambda} : Mortality rates, potentially averaged over
#'     multiple cells in the problem.  Lambda's average over cells according
#'     to the model, M's do not.  M's apply to single cells.
#'
#'     \item coefficients : One or more coefficients in the log model for
#'     lambda.
#'   }
#'
#'   \item \code{intervals} : a matrix containing two-tailed posterior confidence
#'   intervals.  Same dimension as \code{estimates} but columns are lower
#'   and upper endpoints of the intervals. Confidence level of all intervals
#'   is stored in \code{conf.level}.
#'
#'   \item \code{out} : a \code{mcmc.list} object containing the output
#'   of the JAGS runs. This is the raw output from the \code{coda.samples} function
#'   of the \code{coda} package and can be used to check convergence, etc. This
#'   object can be subsetted to chains for one or more parameter by thinking of
#'   the list as a 2-D matrix and using named dimensions.  For example,
#'   \code{x$out[,c("(Intercept)","year")]} produces an mcmc.list object containing only
#'   chains for parameters named "(Intercept)" and "year".  See examples for some
#'   useful plots.
#'
#'   \item \code{jags.model} : a \code{jags.model} object used to estimate the model.
#'   This object can be used to compute additional convergence and model fit statistics,
#'   such as DIC (see \code{\link{rjags::dic.samples}}).
#'   
#'   \item \code{infoCrits} : a named numeric vector containing the WAIC 
#'   and the deviance computed for the model. If \code{computeIC==FALSE},
#'   then \code{infoCrits} is NULL.
#'
#'   \item \code{priors} : the mean and standard deviation of the normal
#'   prior distributions for all coefficients. This differs from the input
#'   \code{priors} parameter because matching of terms and row names has been done.
#'   There is exactly one row in the output \code{priors} for every parameter,
#'   while that is not necessarily true for the input.
#'
#'   \item \code{seeds} : the vector of MCMC random number seeds actually used
#'   in the analysis.  Specify these as input seeds to repeat, exactly, the MCMC
#'   sampling.
#'
#'   \item \code{offset} : the vector of offset terms used in the linear part of the
#'   model.  If no offset terms were specified, this vector contains zeros.
#'
#'   \item \code{coef.labels} : character vector containing the labels of coefficients
#'   in the model.  This is used to distinguish coefficients from derived parameters.
#'   All parameters in \code{ests} not listed here are considered derived. See \code{\link{labels.eoar}},
#'   and \code{\link{coef.eoar}}.
#'
#'   \item \code{call} : the call that invoked this function.
#'
#'   \item \code{data} : the input data set.  This can be handy for prediction and model
#'   selection.
#'
#'   \item \code{converged} : Logical value for whether this routine thinks the
#'   MCMC chain has converged.  The MCMC sampling is deemed to have converged if all Gelman
#'   R-hats are less than some criterion. This function simply calls
#'   \code{\link{checkIsConverged}} with an mcmc.list object containing coefficients only
#'   and a criterion of 1.1.  Convergence checking is only done on lambda model coefficients.
#'
#'   \item \code{Rhats} : Gelman's R-hats for every parameter.  This is one of the output
#'   components of \code{\link{checkIsConverged}}.
#'
#'   \item \code{autoCorrelated} :  Logical value indicating whether this routine
#'   thinks the MCMC is autocorrelated. Normally, one does not want autocorrelation
#'   in the final chains. This function simply calls \code{\link{checkIsAutocorrelated}}
#'   with the mcmc.list and \code{criterion=0.4} and \code{lag=2}.
#'
#'   \item \code{autoCorrs} : a vector of autocorrelations used the judge whether
#'   the chains are autocorrelated.  This is a vector with same length as number of
#'   coefficients, and is the second component of the output of \code{\link{checkIsAutocorrelated}}.
#'   Only model coefficients are checked for autocorrelation.
#'
#'   \item \code{conf.level} : the two-tailed confidence level for all confidence intervals
#'   in \code{intervals}.  Because the entire mcmc.list is returned, one can change
#'   confidence levels without re-running by executing \code{apply(as.matrix(x$out),
#'   2, quantile, p=c(<new quantiles>))}
#'
#' }
#'
#' @author Trent McDonald
#'
#'
#' @seealso \code{\link{labels.eoar}}, \code{\link{coef.eoar}},
#' \code{\link{predict.eoar}}, \code{\link{model.matrix.eoar}}
#'
#' @examples
#' # A 3 year study of 7 sites. 21 "cells". lambda change = 20/year
#' set.seed(9430834) # fixes Y and g of this example, but not the RNG's used in chains
#' ns <- 3
#'
#' ny <- 7
#' g <- data.frame(
#'  alpha = rnorm(ns*ny,70,2),
#'  beta = rnorm(ns*ny,700,25)
#' )
#' Y <- rbinom(ns*ny, c(rep(20,ny), rep(40,ny), rep(60,ny)), g$alpha/(g$alpha+g$beta))
#'
#' df <- data.frame(year=factor(c(rep("2015",ny),rep("2016",ny),rep("2017",ny))),
#'    Year=c(rep(1,ny),rep(2,ny),rep(3,ny)))
#'
#' # Uninformed eoar (use low number of iterations because it's and example)
#' eoar.1 <- eoar(Y~year, g, df, nburn = 1000, niters= 50*10, nthins = 10 )
#'
#' # Repeat and get exact same answers
#' eoar.1 <- eoar(Y~year, g, df, nburn = 1000, niters= 50*10, nthins = 10, seeds=eoar.1$seeds )
#'
#' # Run Informed EoAR.  Assume prior annual lambda estimates
#' # are 10, 15, and 20, all with sd=.5.
#' prior <- data.frame(mean=c(log(10),log(15/10),log(20/10)),
#'   sd=c(0.5,0.5,0.5))
#' row.names(prior) <- c("(Intercept)","year2016","year2017")
#'
#' ieoar <- eoar(Y~year, g, df, priors=prior, nburn = 1000, niters= 50*10, nthins = 10 )
#'
#' # The above chains do not converge due to the low number of iterations used here.
#' # To check convergence and autocorrelation
#'
#' gelman.diag(ieoar$out) # gelmanStats
#' gelman.plot(ieoar$out) # gelmanPlot
#'
#' # Nice traceplots
#' library(lattice)
#' plot(ieoar$out) # tracePlot, all parameters
#' xyplot(ieoar$out[,c("(Intercept)","year2016","year2017")]) # nicer trace of coefficients
#'
#' # Autocorrelation functions
#' acfplot(ieoar$out[,c("(Intercept)","year2016","year2017")], ylim=c(-.2,1), lag.max=300)
#' acfplot(ieoar$out[,c("(Intercept)","year2016","year2017")], ylim=c(-.2,1), thin=2)
#'
#' # Density plots
#' densityplot(ieoar$out[,c("(Intercept)","year2016","year2017")])
#' densityplot(ieoar$out[,c("lambda[1]","lambda[8]","lambda[15]")]) # Mean lambda each year
#'
#' # Correlation among coefficients
#' levelplot(ieoar$out[,c("(Intercept)","year2016","year2017")][[1]])
#'
#' # QQ plots
#' qqmath(ieoar$out[,c("(Intercept)","year2016","year2017")])
#'
#'
eoar <- function(lambda, beta.params, data, offset,
                priors=NULL,
                conf.level=0.9, nburns = 500000, niters = 20000,
                nthins = 10, nchains = 3, nadapt = 3000,
                computeIC = FALSE,
                quiet=FALSE, seeds=NULL, doEoA = FALSE,
                vagueSDMultiplier = 100){
  
  rjags::load.module("dic", quiet = TRUE)

  ## ---- lambdaModel ----
  # Resolve formula for lambda
  if (missing(data))
    data <- environment(lambda)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("lambda", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  names(mf)[names(mf)=="lambda"] <- "formula"
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  Y <- model.response(mf,"any")
  lambda.covars <- if (!is.empty.model(mt)){
    model.matrix(mt, mf, contrasts)
  }
  offset <- as.vector(model.offset(mf))
  ncovars <- ncol(lambda.covars)
  vnames<-dimnames(lambda.covars)[[2]]

  ## ---- initialize ----
  # Make sure one beta dist'n parameter per row
	nyrs <- length(Y)
	if( length(beta.params$alpha) == 1 ){
		alpha.vec <- rep(beta.params$alpha, nyrs)
		beta.vec <- rep(beta.params$beta, nyrs)
	} else {
		alpha.vec <- beta.params$alpha
		beta.vec <- beta.params$beta
	}

	if( length(alpha.vec) != nyrs ) stop("Length of alpha and beta vectors must be 1 or equal length of Y")
	if( length(alpha.vec) != length(beta.vec) ) stop("Lengths of alpha and beta inputs must be equal") # Actually, this can't happen if beta.params is a data frame. Oh well, leave it.


  ## ---- resolvePriors ----
  # Use vague normal priors for coefficients (huge SE's) by default
  sd.n.start <- compVagueSd(Y,alpha.vec,beta.vec,lambda.covars, range.multiplier = vagueSDMultiplier)
  coefTaus <- sd.n.start$vagueSd

  coefMus <- rep(0,ncovars)
  names(coefMus)<-vnames
  if(is.vector(priors)){
    # Use informed prior for intercept parameter
    if(!all(c("mean","sd") %in% names(priors))){
      stop("Priors must contain names 'mean' and 'sd'.")
    } else {
      coefMus["(Intercept)"] <- priors["mean"]
      coefTaus["(Intercept)"] <- priors["sd"]
    }
  } else if(!is.null(priors)){
    # Use informed priors for all coefficience, priors must be a data.frame
    if( !is.data.frame(priors)){
      stop("Priors must be NULL, a vector, or a data.frame.")
    } else if(!all(c("mean","sd") %in% names(priors))){
      stop("Priors data.frame must contain 'mean' and 'sd'.")
    } else if(!any(row.names(priors) %in% vnames)){
      warning(paste("No row.names(priors) matched parameter names. Vague priors used."))
    } else {
      fnd <- names(coefMus) %in% rownames(priors)
      coefMus[fnd] <- priors[names(coefMus)[fnd], "mean"]
      coefTaus[fnd] <- priors[names(coefTaus)[fnd], "sd"]
    }
  }

  # Recall: tau of dnorm in JAGS is 1/variance
  # cat("Prior mean and standard error:\n")
  # cat(paste("mean =", coefMus, "\n"))
  # cat(paste("sd   =", coefTaus, "\n"))

  coefTaus <- 1/coefTaus^2


  ## ---- resolveOffset ----
  if(is.null(offset)){
    offset <- rep(0,nyrs)
  } else {
    if (length(offset) != nyrs)
      stop(gettextf("Length of offset is %d, but should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }


	## ---- bayesModelCode ----
  if(doEoA == FALSE) {
    jagsModel <- "model{

		# Priors
    for(i in 1:ncovars){
      a[i] ~ dnorm( coefMus[i], coefTaus[i] )
    }
    sigmaLambda ~ dunif(0, 100)
    tau <- 1/sigmaLambda^2


    # functional relations
    for(i in 1:nx){
      for(j in 1:ncovars){
        logl[i,j] <- a[j]*lambda.covars[i,j]
      }
      offlink[i] <- exp(offset[i])
      
      lambdaMu[i] <- sum(logl[i,])
    }

		# Likelihood
		for( i in 1:nx ){
			g[i] ~ dbeta(alpha[i], beta[i])
			lambdaLog[i] ~ dnorm(lambdaMu[i], tau)
			lambda[i] <- exp(lambdaLog[i])
			M[i] ~ dpois(offlink[i]*lambda[i])
			Y[i] ~ dbin(g[i], M[i])
		}

		Mtot <- sum(M[])

		}
	"
  } 
  if(doEoA == TRUE) {
    jagsModel <- "model{

		# Priors
    for(i in 1:ncovars){
      a[i] ~ dweib(0.7, 0.03981072)  # informed EoA-style prior
    }


    # functional relations
    for(i in 1:nx){
      for(j in 1:ncovars){
        logl[i,j] <- a[j]*lambda.covars[i,j]
      }
      offlink[i] <- exp(offset[i])
      lambdaMu[i] <- sum(logl[i,])
    }

		# Likelihood
		for( i in 1:nx ){
			g[i] ~ dbeta(alpha[i], beta[i])
			M[i] ~ dpois(offlink[i]*lambdaMu[i])
			Y[i] ~ dbin(g[i], M[i])
		}

		Mtot <- sum(M[])

		}
	"
  }
	

	JAGS.data.0 <- list ( Y = Y,
													nx = nyrs,
													ncovars = ncovars,
													coefTaus = coefTaus,
													coefMus = coefMus,
													alpha = alpha.vec,
													beta = beta.vec,
													offset = offset,
													lambda.covars = lambda.covars)


  writeLines(jagsModel, "model.txt")
  #cat(jagsModel)


  ## ---- initialValues ----
  # MCMC sample size settings:


  if( is.null(seeds) ){
    tmp.mult <- 10^(8)
    seeds <- round(runif(nchains,0,tmp.mult))
  } else if (length(seeds) != nchains) {
    stop(paste("MCMC random number seeds should be NULL or length", nchains ))
  }

  Inits <- function(x,strt,seed){
  		gg <-rbeta(x$nx, x$alpha, x$beta)
  		M <- ceiling(x$Y / gg) + 1
  		ys <- unique(x$Y)
  		a <- rep(0,x$ncovars)
  		a[1] <- log(mean(M)) + 1E-10  # to avoid true zero and errors
  		list ( a = a,
  					 M = M,
  					 g=gg,
  					 .RNG.name="base::Mersenne-Twister",
  					 .RNG.seed=seed )
  }

  inits<-vector("list",nchains)
  for(i in 1:nchains){
    inits[[i]]<-Inits(JAGS.data.0, sd.n.start, seeds[i])
  }


  # Parameters to be monitored by WinBUGS
  params <- c("a", "M", "lambda", "Mtot", "sigmaLambda")


  ## ---- jagsRun ----

  # Initialize the chains and adapt
  (t1=Sys.time())

  assign("JAGS.data.0", JAGS.data.0, envir=.GlobalEnv )
  assign("inits", inits, envir=.GlobalEnv )

  jags = jags.model(file="model.txt",
  									data=JAGS.data.0,
  									inits=inits,
  									n.chains=nchains,
  									n.adapt=nadapt,
  									quiet = quiet)

  #   Run out in the chain a ways
  if(!quiet) cat("Burnin...\n")
  out = update(jags, n.iter=nburns, progress.bar=ifelse(quiet, "none","text"))

  # Run the MCMC chains
  if(!quiet) cat("MCMC Sampling...\n")
  out = coda.samples(jags,
  									 variable.names=params,
  									 n.iter=niters,
  									 thin=nthins,
  									 progress.bar=ifelse(quiet, "none","text"))
  
  ics <- NULL
  if(computeIC) {
    if(!quiet) cat("Computing WAIC and deviance...\n")
    s <- jags.samples(jags, 
                      c("deviance", "WAIC"), 
                      type = "mean", 
                      n.iter = niters,
                      thin = nthins,
                      progress.bar=ifelse(quiet, "none","text")) 
    s <- lapply(s, unclass)
    ics <- sapply(s, sum)
  } 
   
  (t2=Sys.time())
  t3 <- t2-t1
  if(!quiet) cat(paste("Execution time:", round(t3,2), attr(t3,"units"), "\n\n"))

  ## ---- mcmcChecking ----
  #print(summary(out))

  # extract the coefficients and label em
  # Careful.  You can't have any other parameters that start with 'a'.
  out.coefs <- out[,grep("^a",varnames(out))]
  varnames(out)[grep("^a",varnames(out))] <- vnames
  varnames(out.coefs) <- vnames

  # Check convergence
  conv <- checkIsConverged(out.coefs,1.1,quiet)

  # Check autocorrelation
  auto <- checkIsAutocorrelated(out.coefs, 0.4, 2, quiet)



  ## ---- quantiles ----
  alpha <- c(0,0.5,1) + c(1,0,-1)*(1-conf.level)/2 # median plus two extreme-er quantiles
  out.sumry <- summary(out, quantiles=alpha)
  out.sd <- out.sumry$statistics[,c("SD")]
  out.sumry <- out.sumry$quantiles
  out.sd <- cbind(out.sumry[,"50%"], out.sd)
  dimnames(out.sd)[[2]] <- c("Estimate","SD")
  out.sumry <- list(estimates=out.sd,
                    intervals=out.sumry[,-2]
                    )

  ## ---- Done ----
  priors.df <- data.frame(mean=coefMus, sd=1/sqrt(coefTaus))
  ans <- c(
    out.sumry,
    list(
    out=out,
    jags.model=jags,
    infoCrits = ics,
    priors = priors.df,
    seeds=seeds,
    offset=offset,
    coef.labels=vnames,
    call=cl,
    data=data
    ),
    conv,
    auto,
    conf.level=conf.level,
    terms = mt
  )

  class(ans) <- "eoar"
  ans

}



