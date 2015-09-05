\name{spm-package}
\alias{spm-package}
\alias{spm}
\docType{package}
\title{
Stochastic Process Modeling
}
\description{
Stochastic Process Modeling (SPM) perform estimation of parameters of some stochastic process of interest. This package is mainly for estimating parameters from clinical traits.
}
\details{
\tabular{ll}{
Package: \tab spm\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-09-05\cr
License: \tab GPL\cr
}

}
\author{
Ilya Y. Zhbannikov

Maintainer: Ilya Y. Zhbannikov <i.zhbannikov@mail.ru>
}
\references{

Anatoli I. Yashin, Konstantin G. Arbeev, Igor Akushevich, Aliaksandr Kulminski, Lucy Akushevich, Svetlana V. Ukraintseva, Stochastic model for analysis of longitudinal data on aging and mortality, Mathematical Biosciences, Volume 208, Issue 2, August 2007, Pages 538-551, ISSN 0025-5564, http://dx.doi.org/10.1016/j.mbs.2006.11.006.

Anatoli I. Yashin , Konstantin G. Arbeev, Aliaksandr Kulminski, Igor Akushevich, Lucy Akushevich, Svetlana V. Ukraintseva, Health decline, aging and mortality: how are they related? Biogerontology
June 2007, Volume 8, Issue 3, pp 291-302
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package
SPM
Stochastic
longevity
clinical
medicine }
\examples{
library(spm)
# Reading longitude data:
longdat <- read.csv(system.file("data","longdat.csv",package="spm"))
# Prepare data for optimization:
vitstat <- read.csv(system.file("data","vitstat.csv",package="spm"))
# Remove unneeded NAs:
longdat.nonan <- longdat[which(is.na(longdat$Age) == F),]
vitstat.nonan <- vitstat[which(is.na(vitstat$BirthCohort) == F),]
data=prepare_data(longdat=longdat.nonan, vitstat=vitstat.nonan,interval=1, col.status="IsDead", col.id="ID", col.age="Age", col.age.next="AgeNext", col.age.event="LSmort", covariates=c("DBP"), verbose=T)
# Parameters estimation:
pars=spm(data,k = 1)
pars
}
