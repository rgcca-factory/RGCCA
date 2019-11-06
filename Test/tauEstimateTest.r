# Tau estimate

x=matrix(rnorm(12),3,4)
tau.estimate(x,na.rm=TRUE)
tau.estimate(x,na.rm=FALSE)

y=x
y[1,2]=NA
tau.estimate(y,na.rm=TRUE)
