# Evidence of Absence Regression Package

These are routines implement Evidence of Absence Regression (EoAR). 
**EoAR** relates the number of found entities after a series of
searches to covariates to estimate the number of missed entities. 
Special cases are 
the **Evidence of Absence (EoA)** model of Huso et al. (2015) and the
**Informed Evidence of Absence (IEoA)** approaches.

## How to git it:

To obtain a local copy of the source, download it from GitHub using the following 
commands in a git shell: 
```
cd <directory you want>  
git clone https://github.com/tmcd82070/eoar.git  
```

The above commands will download all source from GitHub to your computer.  
Among other things, 
you should see a `DESCRIPTION` file and `R` directory.  

**You do not need a local copy** unless you plan to edit the code.  If you 
simply want to run `EoAR`, use the 'devtools and GitHub' installation method below. 


## Intalling:

`EoAR` depends on `rjags`.  It is best to install `rjags` prior to 
atempting install of `EoAR`.  Prior to attempting installation of 
`rjags`, it is best to install JAGS.  If you already have the 
JAGS runtime on your machine, skip this step and install 
`rjags` immediately. 

### Install JAGS run time

To install the JAGS 
runtime, navigate here: 
`http://www.sourceforge.net/projects/mcmc-jags/files`, and 
download the JAGS installer. Execute it and accept all defaults. 

### Install `rjags`
Issue the following command in R:
```
install.packages("rjags")
```

You have JAGS and `rjags` installed properly if you attach `rjags` and 
do not get an error message.  I.e., something like the following:
```
> library(rjags)
Loading required package: coda
Linked to JAGS 4.3.0
Loaded modules: basemod,bugs
```

### Install `EoAR`

#### Using `devtools` and GitHub

You can install the latest release of `EoAR` directly from GitHub without 
obtaining a local copy.  Issue the following command:
```
devtools::install_github("tmcd82070/EoAR")
```

#### Using `devtools` and local repository

If you cloned the `EoAR` repository, you have a local copy of the source on your
hard drive.  Open R and `setwd()` to the directory containing the `DESCRIPTION` file. In R issue the following:
```
library(devtools)  
document()  
install()   
```

#### Manual install

Open a command window, change directory to the folder containing `DESCRIPTION` and issue 
the following command: 
```
r CMD INSTALL EoAR
```

## To Contribute

If you change something, and it's useful, I would be very interested to hear about it. 

## Usage Example

The main routine is `eoar`.  It takes a count vector, model for lambda, and g-values, 
Here is an example : 

This is fake data from a three year study on seven sites.  First, the 
alpha and beta parameters for g-value distributions, one per year.   
```
ns <- 3  
ny <- 7  
g <- data.frame(  
  alpha = rnorm(ns*ny,70,2),  
  beta = rnorm(ns*ny,700,25)  
)
```

The following command generates a carcasses count vector, one count per site per year.  
```
Y <- rbinom(ns*ny, c(rep(20,ny), rep(40,ny), rep(60,ny)), g$alpha/(g$alpha+g$beta))
```

The following command generates a fake covariate data frame.  
This example data frame contains *Year* as a linear 
effect (1,2,3,etc) and *year* as a factor (2015, 2016, 2017, etc).  

```
df <- data.frame(year=factor(c(rep("2015",ny),rep("2016",ny),rep("2017",ny))),  
    Year=c(rep(1,ny),rep(2,ny),rep(3,ny)))
```

The following relates carcass deposition rates to *year* using 
vague priors for coefficients:     
```
eoa.1 <- eoar(Y~year, g, df )
```
The following uses informed distributions:

```
# Assume prior mean is 10 and prior sd is 3  
# Fit intercept-only model to get one mean lambda   
intMean <- 2*log(10) - 0.5*log(3^2 + 10^2)  
intSd <- sqrt(-2*log(10) + log(3^2 + 10^2))  
prior <- data.frame(mean=intMean, sd=intSd, row.names="(Intercept)")  
eoa.1 <- eoa(Y~1, g, df, priors=prior )  
```

After either run, you should check convergence.  
To do so, run a traceplot and Gelman stats.  Any Rhats > 1.1 indicate suspect 
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
