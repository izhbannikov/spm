% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spm.impute.R
\name{spm.impute}
\alias{spm.impute}
\title{Multiple Data Imputation with SPM}
\usage{
spm.impute(x, id = 1, case = 2, t1 = 3, t2 = 3, covariates = 4,
  minp = 5, theta_range = seq(0.01, 0.2, by = 0.001))
}
\arguments{
\item{x}{A longitudinal dataset with missing observations}

\item{id}{A name (text) or index (numeric) of ID column. Default: 1}

\item{case}{A case status column name (text) or index (numeric). Default: 2}

\item{t1}{A t1 (or t if short format is used) column name (text) or index (numeric). Default: 3}

\item{t2}{A t2 column name (if long format is used) (text) or index (numeric). Default: 4}

\item{covariates}{A list of covariate column names or indices. Default: 5}

\item{minp}{Number of imputations. Default: 5}

\item{theta_range}{A range of parameter theta used for optimization, default: seq(0.01, 0.15, by=0.001).}
}
\value{
A list(imputed, imputations)

imputed An imputed dataset.

imputations Temporary imputed datasets used in multiple imputaitons.
}
\description{
Multiple Data Imputation with SPM
}
\examples{
\dontrun{
library(stpm) 
##Data preparation ##
data <- simdata_discr(N=1000, dt = 2)
miss.id <- sample(x=dim(data)[1], size=round(dim(data)[1]/4)) # ~25\% missing data
incomplete.data <- data
incomplete.data[miss.id,5] <- NA
incomplete.data[miss.id-1,6] <- NA
## End of data preparation ##

# Estimate parameters from the complete dataset #
p <- spm_discrete(data, theta_range = seq(0.075, 0.09, by=0.001))
p

##### Multiple imputation with SPM #####
imp.data <- spm.impute(x=incomplete.data, 
                      minp=5, 
                      theta_range=seq(0.075, 0.09, by=0.001))$imputed
head(imp.data)
## Estimate SPM parameters from imputed data and compare them to the p ##
pp.test <- spm_discrete(imp.data, theta_range = seq(0.075, 0.09, by=0.001))
pp.test
}
}
