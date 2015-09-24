# Stochastic Process Modeling (SPM)
## Features
* Data simulation
** Continuous
** Discrete
* Optimisation
** Continuous
** Discrete
** Time-dependant coefficients
## How to use
```
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
```