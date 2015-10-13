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
data=prepare_data(longdat=longdat.nonan, vitstat=vitstat.nonan,interval=1, col.status="IsDead", col.id="ID", col.age="Age", col.age.event="LSmort", covariates=c("DBP"), verbose=T)
# Parameters estimation:
pars=spm(data,k = 1)
pars
```