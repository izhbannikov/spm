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

**Note:** for Windows, please install Rtools: https://cran.r-project.org/bin/windows/Rtools/
**Note:** compilation under Windows may fail if you run Windows on a virtual machine!
Then:
```
install.packages("devtools")
library(devtools)
install_github("izhbannikov/spm")
```

## Examples

### Test discrete simulation
```
library(stpm)
## Test discrete
data <- simdata_discr(N=1000)
pars <- spm_discrete(data)
pars
```

### Test projection

```
library(stpm)
# Setting up the model
model.par <- list()
model.par$a <- matrix(c(-0.05, 1e-3, 2e-3, -0.05), nrow=2, ncol=2, byrow=TRUE)
model.par$f1 <- matrix(c(90, 35), nrow=1, ncol=2)
model.par$Q <- matrix(c(1e-8, 1e-9, 1e-9, 1e-8), nrow=2, ncol=2, byrow=TRUE)
model.par$f <- matrix(c(80, 27), nrow=1, ncol=2)
model.par$b <- matrix(c(6, 2), nrow=2, ncol=2)
model.par$mu0 <- 1e-6
model.par$theta <- 0.09
# Projection
# Discrete-time model
data.proj.discrete <- spm_projection(model.par, N=5000, ystart=c(80, 27), tstart=c(30, 60))
plot(data.proj.discrete$stat$srv.prob)
# Continuous-time model
data.proj.continuous <- spm_projection(model.par, N=5000, ystart=c(80, 27), model="continuous", gomp=TRUE)
plot(data.proj.continuous$stat$srv.prob)
```

### Time-dependent model
```
library(stpm)
model.par <- list(at="-0.05", f1t="80", Qt="2e-3", ft="80", bt="15", mu0t="1e-5*exp(0.08*t)")
data.proj.time_dependent <- spm_projection(model.par, N=500, ystart=80, model="time-dependent", sd=4)
plot(data.proj.time_dependent$stat$srv.prob)
```

### Test prepare_data()
```
library(stpm)
data <- prepare_data(x=system.file("data","longdat.csv",package="stpm"), y=system.file("data","vitstat.csv",package="stpm"))
head(data[[1]])
head(data[[2]])
```

### Test simdata_gen_cont()

```
library(stpm)
dat <- simdata_gen(N=500)
head(dat)
```

### Test spm_gen()
```
library(stpm)
#Reading the data:
data <- simdata_gen(N=100)
head(data)
#Parameters estimation:
pars <- spm_gen(gendat=data)
pars
```

### Test simdata_cont()
```
library(stpm)
dat <- simdata_cont(N=50)
head(dat)
```

### Test simdata_discr()
```
library(stpm)
data <- simdata_discr(N=100)
head(data)
```

### Test simdata_time-dep()
```
library(stpm)
dat <- simdata_time_dep(N=100)
head(dat)
```

### Test spm_continuous()
```
library(stpm)
data <- simdata_cont(N=50)
pars <- spm_continuous(dat=data,a=-0.05, f1=80, Q=2e-8, f=80, b=5, mu0=2e-5)
pars
```

### Test spm_discrete()
```
library(stpm)
data <- simdata_discr(N=10)
pars <- spm_discrete(data)
pars
```

### Test spm()
```
library(stpm)
data.continuous <- simdata_cont(N=1000)
data.discrete <- simdata_discr(N=1000)
data <- list(data.continuous, data.discrete)
p.discr.model <- spm(data)
p.discr.model
p.cont.model <- spm(data, model="continuous")
p.cont.model
p.td.model <- spm(data, model="time-dependent",formulas=list(at="aa*t+bb", f1t="f1", Qt="Q", ft="f", bt="b", mu0t="mu0"), start=list(a=-0.001, bb=0.05, f1=80, Q=2e-8, f=80, b=5, mu0=1e-3))
p.td.model
```
