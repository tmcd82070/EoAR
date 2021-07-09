# Evidence of Absence Regression Package

These routines implement Evidence of Absence Regression (EoAR) methods 
described in McDonald et al. (2021, *Evidence of absence regression: a binomial N-mixture model for estimating  
fatalities at wind energy facilities*, Ecological Applications, In press, URL pending). 
The **EoAR** method estimates the number of (found + missed) entities after a 
series of searches by using probability of detection and 
covariate relationships. 
Special cases are 
the **Evidence of Absence (EoA)** model of Huso et al. (2015) and the
**Informed Evidence of Absence (IEoA)** approaches.

# Data

Example data analyzed in McDonald et al., (2021) are available on the 
Dryad digital repository: [https://doi.org/10.5061/dryad.2rbnzs7jh](https://doi.org/10.5061/dryad.2rbnzs7jh).

# Installation

**EoAR** uses MCMC estimation methods implemented in JAGS.  Hence, **EoAR** requires 
two prerequisites. 

## 1. Install the JAGS run time

The JAGS runtime is the program that actually performs the MCMC computations. 
Install the JAGS runtime by navigating here: [http://www.sourceforge.net/projects/mcmc-jags/files](http://www.sourceforge.net/projects/mcmc-jags/files).
Download the JAGS installer. Execute it and accept all defaults. 

## 2. Install the `rjags` R package

The `rjags` R package communicates with the JAGS runtime and makes running MCMC code in 
R convienient. `rjags` is available on R's CRAN repository. To install `rjags` issue the following command in R:
```
install.packages("rjags")
```

**Test the JAGS installation**: You can test the JAGS and `rjags` installation by simply attaching 
`rjags` in R. If JAGS and `rjags` are okay, attaching `rjags` using the `library` 
command will result in output similar to the following:
```
> library(rjags)
Loading required package: coda
Linked to JAGS 4.3.0
Loaded modules: basemod,bugs
```

If `rjags` cannot find the JAGS runtime, an error will be issued after the 
`library(rjags)` command. 

## 3. Install the `EoAR` package

The easist way to install `EoAR` is directly from GitHub using 
routines in the `devtools` package.  If you do not have `devtools`, install 
it with the following command: 
```
install.packages("devtools")
```
Then, to install `EoAR`, issue the following command:

```
devtools::install_github("tmcd82070/EoAR")
```

# Usage Example

The main routine is `eoar`.  The `eoar` routine takes a count vector, a model for lambda, and g-values 
as inputs. Following is an `eoar` example on fake (simulated) data. 

The following generates fake data for a three year study on seven sites.  We first generate 
*alpha* and *beta* parameters from common normal distributions which we later 
use to generate site and year specific g-values.
```
ns <- 7  # Number of sites
ny <- 3  # Number of years
g <- data.frame(  
  alpha = rnorm(ns*ny,70,2),  
  beta = rnorm(ns*ny,700,25)  
)
```

In this example, we assume that the true average number of carcasses per site 
increases each year, but does not vary by site.  We assume the true average number 
of carcasses per site in year 1 is 20, year 2 is 40, and year 3 is 60. 
We assume carcasses at each site each year are detected with probabilities
equal to $\alpha / (\alpha + \beta)$, which are the means of beta distributions 
assumed for g-values. 
The following generates random observed carcasses counts, one per site per year.  
```
meanYr1 <- 20
meanYr2 <- 40
meanYr3 <- 60
Y <- rbinom(ns*ny, c(rep(meanYr1,ns), rep(meanYr2,ns), rep(meanYr3,ns)), g$alpha/(g$alpha+g$beta))
```

In this example, we fit a linear trend and annual categories to the true number of carcasses.  
We construct the linear *Year* and factor *year* covariates in a data frame using the 
following code:  

```
df <- data.frame(year=factor(c(rep("2015",ns),rep("2016",ns),rep("2017",ns))),  
    Year=c(rep(1,ns),rep(2,ns),rep(3,ns)))
```

Finally, we relate true carcass deposition rates to *year* and *Year*. In this 
simple example, we assume that we have correctly estimated all 
g-values by using the generated $\alpha$ and $\beta$ parameters.  
This *EoAR* run uses vague priors for the coefficients.     
```
eoa.1 <- eoar(Y~year, g, df )
eoa.2 <- eoar(Y~Year, g, df )
```

## Informed Priors

When appropriate, it is possible to inform the *EoAR* coefficient's
prior distributios. 
The following assumes that the prior mean number of carcasses per 
site is 10 with a standard deviation of 3. The following code 
fits an intercept-only model.

```
intMean <- 2*log(10) - 0.5*log(3^2 + 10^2)  
intSd <- sqrt(-2*log(10) + log(3^2 + 10^2))  
prior <- data.frame(mean=intMean, sd=intSd, row.names="(Intercept)")  
eoa.1 <- eoa(Y~1, g, df, priors=prior )  
```

## Model Checks

After running *EoAR*, you should check convergence.  
We suggest running traceplots and Gelman stats.  Any Rhats > 1.1 indicate suspect 
convergence. The following commands are useful for inspecting 
mixing and convergence:
```
library(lattice)
xyplot(ieoa.1$out[,labels(ieoa.1)])
acfplot(ieoa.1$out[,labels(ieoa.1)])   
densityplot(ieoa.1$out[,labels(ieoa.1)])  
gelman.diag(ieoa.1$out) # gelmanStats  
gelman.plot(ieoa.1$out) # gelmanPlot  
```
