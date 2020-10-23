
rm(list = ls()) 
library(quantregForest)
library(stringr)
library(devEMF)
File="DGP1 (Sample)"
set.seed(1)

# Define a function ci.qcm of T0 to compute the coverage rate for T1 posttreatment periods

# Other parameters and DGP are all fixed 

# ci.qcm refers to "confidence interval for quantile control method"

trTime=21
trId=1
T=21
F_K=3
C_K=2
G_K=3
K=2
N=10
X=array(dim=c(N,T,K))
F=array(dim=c(T,F_K))
C=array(dim=c(N,C_K))
e=array(dim=c(N,T,K))
G=array(dim=c(N,G_K))
R=array(dim=c(N,K))
Y=array(dim=c(N,T))
u=array(dim=c(N,T))
beta=c(1,2)
F[,1:F_K]=rnorm(n=T*F_K,mean=0,sd=1)
C[,1:C_K]=runif(n=N*C_K,min=1,max=2)
e[,,]=rchisq(n=N*T*K,df=1)-1
G[,]=rnorm(n=N*G_K,mean=0,sd=1)
R[,]=runif(n=N*K,min=0.1,max=0.9)
X[,1,]=runif(n=N*K)
u[,]=rnorm(n=N*T)
for(t in 2:T){
  for(i in 1:N){
    for(k in 1:K){
      X[i,t,k]=1+R[i,k]*X[i,t-1,k]+sum(C[i,]*c(G[i,k],F[t,k]))+e[i,t,k]
    }    
  }  
}

for(i in 1:N){
  for(t in 1:T){
    Y[i,t]=sum(X[i,t,]*beta)+sum(G[i,]*F[t,])+u[i,t]
  }
}
train=1:(trTime-1)
test=trTime:T
Y[trId,test]=Y[trId,test]+1
y=matrix(t(Y[trId,]),ncol=1)
y_train=y[train,]
y_test=y[test,]
x=NULL
for(k in 1:K){
  x=cbind(x,t(X[,,k]))
}
x=cbind(x,t(Y[-trId,]))
x_train=x[train,]
x_test=x[test,]

qrf <- quantregForest(x=x_train, y=y_train,ntree=5000,nodesize=5,replace=FALSE)
conditionalQuantiles=predict(qrf, x_test, what=c(0.025,0.975))
qr025.pred=conditionalQuantiles[,1]
qr975.pred=conditionalQuantiles[,2]

treat_upr=y_test-qr025.pred
treat_lwr=y_test-qr975.pred

if (treat_upr<1 || treat_lwr>1) {
  cover_forest=0
} else{
  cover_forest=1
}

# Median Treatment Effect
conditionalQuantiles=predict(qrf,x,what=c(0.025,0.5,0.975))
qr025.pred=conditionalQuantiles[,1]
qr500.pred=conditionalQuantiles[,2]
qr975.pred=conditionalQuantiles[,3]
tr.median=y-qr500.pred

# Plot
qrf.plot=function(){
  plot(c(1,T),c(-10,30),type="n",xlab="Time",ylab="Value")
  lines(1:T,y,col="blue")
  lines(1:T,qr025.pred,col="red",lty=2)
  lines(1:T,qr500.pred,col="red",lty=1)
  lines(1:T,qr975.pred,col="red",lty=2)
  legend("topright",legend=c("Y","Median.pred","95%CI.pred"),lty=c(2,1,2),col=c("blue","red","red"))
}
qrf.plot()
emf(file=paste(File,"Photo.emf"), width=8, height=6,bg="transparent")
qrf.plot()
dev.off()

# End