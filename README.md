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
#Prepare data for optimization
data <- prepare_data(x=system.file("data","longdat.csv",package="spm"), y=system.file("data","vitstat.csv",package="spm"))
#Parameters estimation (default model: discrete-time):
p.discr.model <- spm(data)
p.discr.model
# Continuous-time model:
p.cont.model <- spm(data, model="continuous")
p.cont.model
# Model with time-dependent coefficients:
data <- prepare_data(x=system.file("data","longdat.csv",package="spm"), y=system.file("data","vitstat.csv",package="spm"), covariates="BMI")
p.td.model <- spm(data, model="time-dependent")
p.td.model
```