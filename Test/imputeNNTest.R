# '# Test imputeNN
# 
# '''
set.seed(42);X1=matrix(rnorm(500),100,5);
set.seed(22);X2=matrix(rnorm(400),100,4);
set.seed(2);X3=matrix(rnorm(700),100,7);
# usual test 
X1[1,]=NA
X2[7,1]=NA
X2[5,1]=NA
X3[3,]=NA
X3[4,]=NA
rownames(X1)=rownames(X2)=rownames(X3)=paste("S",1:100,sep="")
colnames(X1)=paste("A",1:5,sep="")
colnames(X2)=paste("A",6:9,sep="")
colnames(X3)=paste("A",10:16,sep="")
A=list(X1,X2,X3)  
Ares=imputeNN(A,k=1,output="mean", klim=NULL,scale=TRUE,sameBlockWeight=TRUE)
 

# x=1:10
# taille=scale(1:10+rnorm(10,0,0.01))
# taille2=scale(1:10+rnorm(10,0,1))
# poids=scale(c(1,1,1,1,2,2,2,1,2,2))
# lm(x~taille+taille2+poids)
# summary(lm(x~taille+taille2+poids))
# cor(taille,taille2)
# 
# cor(x,taille)
# cor(x,taille2)
# cor(x,poids)