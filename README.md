# How to use
```
library(mice)
library(spm)
# Reading longitude data:
longdat <- read.csv(system.file("data","longdat.csv",package="spm"))
# Remove unneeded NAs:
longdat.nonan <- longdat[which(is.na(longdat$Age) == F),]
# Prepare data for optimization:
ans=prepare_data(longdat.nonan,covariates=c("DBP","BMI"))
# Look at the data
head(ans)
# Analysis
res=spm(ans)
# Look at the results:
res
```