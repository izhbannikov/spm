# Stochastic Process Model (SPM)
## Features
### Data simulation
* Continuous (one- and multiple-dimensions)
* Discrete (one- and multiple-dimensions)

### Optimisation
* Continuous (one- and multiple-dimensions)
* Discrete (one- and multiple-dimensions)
* Time-dependant coefficients (one-dimensional optimisation)

### How to install
#### Stable version
```
install.packages("stpm")
```
#### The most recent version
```
install.packages("devtools")
library(devtools)
install_github("izhbannikov/spm")
```

## How to use
```
library(stpm)
#Prepare data for optimization
data <- prepare_data(x=system.file("data","longdat.csv",package="stpm"), 
				   y=system.file("data","vitstat.csv",package="stpm"))
# Warining! If you experience an error telling that data must be in 
# csv or sas7bdat format,
# use the following commands to prepare the test data:
#install.packages("R.utils")
#library(R.utils)
#gunzip(system.file("data","longdat.csv.gz",package="stpm"))
#gunzip(system.file("data","vitstat.csv.gz",package="stpm"))
#data <- prepare_data(x=system.file("data","longdat.csv",package="stpm"), 
#				   y=system.file("data","vitstat.csv",package="stpm"))
#
#Parameters estimation (default model: discrete-time):
p.discr.model <- spm(data)
p.discr.model
```