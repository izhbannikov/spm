# Stochastic Process Modeling (SPM)
## Features
### Data simulation
* Continuous (one- and multiple-dimensions)
* Discrete (one- and multiple-dimensions)

### Optimisation
* Continuous (one- and multiple-dimensions)
* Discrete (one- and multiple-dimensions)
* Time-dependant coefficients (one-dimensional optimisation)

## How to use
```
library(spm)
# Reading data
## Longitudinal studies:
longdat <- read.csv(system.file("data","longdat.csv",package="spm"))
## Vital statistics:
vitstat <- read.csv(system.file("data","vitstat.csv",package="spm"))
# Prepare data for optimization:
data=prepare_data(longdat=longdat, vitstat=vitstat,interval=1, col.status="IsDead", col.id="ID", col.age="Age", col.age.event="LSmort", covariates=c("DBP"), verbose=T)
# Parameters estimation:
pars=spm(data,k = 1)
pars
```

## SPM for time-dependent coefficients:
```
library(spm)
# Reading data
## Longitudinal studies:
longdat <- read.csv(system.file("data","longdat.csv",package="spm"))
## Vital statistics:
vitstat <- read.csv(system.file("data","vitstat.csv",package="spm"))
# Prepare data for optimization:
data=prepare_data(longdat=longdat, vitstat=vitstat,interval=1, col.status="IsDead", col.id="ID", col.age="Age", col.age.eve$
# Parameters estimation:
ans <- spm_time_dep(data[[1]][,2:6], formulas=list(at="a_y+b_y*t", 
                                                     f1t="a_f1+b_f1*t", 
                                                     Qt="a_q + b_q*t", 
                                                     ft="a_f+b_f*t", 
                                                     bt="b", 
                                                     mu0t="mu_0a*exp(mu_0b*t)"), 
                      start=list(a_y=-0.1, b_y=0.001,
                                 a_f1=80, b_f1=0.002,
                                 a_q=1e-3, b_q=1e-4,
                                 a_f=80, b_f=0.009,
                                 mu_0a=1e-4, mu_0b=0.08,
                                 b=5), 
                      verbose = T)
ans
```
