install.packages("mvtnorm")
library(mvtnorm)
#DGP
M=300
n=2000
gamma.sim=c(0.01,0.01); k=length(gamma.sim)
schaetzung=array(0, c(M,3,k))               
D=array(0, c(k, k, n))
dcc.schaetz=matrix(0,ncol=2,nrow=M)
alpha.sim=c(0.05,0.05)
beta.sim=c(0.9,0.9)
A = 0.05
B = 0.94
sd=matrix(0,ncol=2,nrow=M)
abl1 = matrix(0,ncol=1,nrow=M);abl2 = matrix(0,ncol=1,nrow=M)
var=array(0, c(2, 2, M))
S_12=numeric(M)
S_true = matrix(c(1,0.6,0.6,1),ncol=2,nrow=2)
R_m=numeric(M)


#Funktion fuer DGP mit DCC-Prozess
DGP_funktion = function(n){
  r=matrix(0, ncol=k, nrow=n)
  sigma.sim = array(0, c(k, k, n))
  eps = matrix(0,ncol=n,nrow=k)
  Q = array(0, c(k, k, n))
  R = array(0, c(k, k, n))
  D = array(0, c(k, k, n)); Dinv = array(0, c(k, k, n))
  D[,,1] = matrix(c(sqrt(gamma.sim[1]/(1-alpha.sim[1]-beta.sim[1])),0,0,sqrt(gamma.sim[2]/(1-alpha.sim[2]-beta.sim[2]))),ncol=2,nrow=2) 
  Dinv[,,1] = solve(D[,,1])                
  Q[,,1] = S_true  
  #Datengenerierung
  sigma.sim[,,1] = matrix(c(gamma.sim[1]/(1-alpha.sim[1]-beta.sim[1]),0,0,gamma.sim[2]/(1-alpha.sim[2]-beta.sim[2])),ncol=2,nrow=2)
  for (t in 2:n){
    sigma.sim[,,t] = diag(c(gamma.sim[1],gamma.sim[2])) + diag(c(alpha.sim[1],alpha.sim[2]))*r[t-1,]%*%t(r[t-1,]) + diag(c(beta.sim[1],beta.sim[2]))*sigma.sim[,,t-1] 
    D[,,t] = diag(sqrt(c(sigma.sim[1,1,t], sigma.sim[2,2,t])))
    Dinv[,,t] = solve(D[,,t]) 
    eps[,t-1] = Dinv[,,t-1]%*%r[t-1,]
    Q[,,t] = S_true*(1-A-B) + A*eps[,t-1]%*%t(eps[,t-1]) + B*Q[,,t-1]
    Qast <- diag(c((diag(Q[, , t])^-0.5)),2,2)
    R[, , t] <- Qast %*% Q[, , t] %*% Qast
    r[t,] = rmvnorm(n=1, mean=rep(0,k), sigma=D[,,t]%*%R[,,t]%*%D[,,t]) 
  }
  
  result=list(r,R,Q,t(eps)) #1:r#2:R#3:Q#4:eps
  return(result)
}


#EWMA Modell
H=array(0,c(k,k,n))
h=matrix(0,ncol=k,nrow=n)
R_bias=numeric(M)
R_12=numeric(n)
for (i in 1:M){

DGP=DGP_funktion(n)
r=matrix(c(unlist(DGP[1])),ncol=k,nrow=n)
R_true=array(c(unlist(DGP[2])),c(k,k,n))
R_true=R_true[1,2,]
lambda=0.94
#Initialwerte
H[,,1]=cov(r)                           #Kovarianzmatrix
h[1,1]=var(r[,1]);h[1,2]=var(r[,2])     #Varianzen
R_12[1]=H[1,2,1]/(sqrt(h[1,1])*sqrt(h[1,2]))
for (t in 2:n){
  for(u in 1:k){
h[t,u]=lambda*h[t-1,u] + (1-lambda)*r[t-1,u]^2 
}
H[,,t]=(1-lambda)*(r[t-1,]%*%t(r[t-1,])) + lambda*H[,,t-1]
R_12[t]=H[1,2,t]/(sqrt(h[t,1])*sqrt(h[t,2]))
}
R_bias[i]=mean(abs(R_12-R_true))
}
R_MAE=mean(R_bias)
#Output
R_MAE