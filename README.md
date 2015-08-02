# How to use
```
library(mice)
library(spm)
ndim = 2 # 2-dimensional optimization

# Reading longitude data:
longdat <- read.csv(system.file("data","longdat.csv",package="spm"))


# Prepare data for optimization:
vitstat <- read.csv(system.file("data","vitstat.csv",package="spm"))
# Remove unneeded NAs:
longdat.nonan <- longdat[which(is.na(longdat$Age) == F),]
vitstat.nonan <- vitstat[which(is.na(vitstat$BirthCohort) == F),]
ans=prepare_data(longdat=longdat.nonan, vitstat=vitstat.nonan,covariates=c("DBP","BMI"))
# Look at the data
head(ans[[1]])
head(ans[[2]])

# Optimization and parameters estimation:
pars=spm(ans[[2]],k = ndim)
# Look at the results:
pars


```