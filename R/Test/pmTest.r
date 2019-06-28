'# Test pm 

'''
wd="/home/caroline.peltier/Bureau/RGCCAtoPull"
setwd(wd)
source("pm.r")
set.seed(42); A=matrix(rnorm(15),3,5)
set.seed(34); B=matrix(rnorm(20),5,4)

# loading the new function output
T1=Sys.time();pmToGet=A%*%B;T2=Sys.time(); Tdiff1=T2-T1
T1=Sys.time();pmRes=pm(A,B);T2=Sys.time();Tdiff2=T2-T1
# testing the same results than previously
all.equal(pmToGet,pmRes)
Tdiff2-Tdiff1       
# loading the new function output
A[3,5]=0
B[1,1]=0
T1=Sys.time();pmToGet=A%*%B;T2=Sys.time(); Tdiff1=T2-T1
A[3,5]=NA
B[1,1]=NA
T1=Sys.time();pmRes=pm(A,B);T2=Sys.time();Tdiff2=T2-T1
all.equal(pmToGet,pmRes)
Tdiff2-Tdiff1

