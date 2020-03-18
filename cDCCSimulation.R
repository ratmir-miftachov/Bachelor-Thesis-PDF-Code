#Simulation cDCC
install.packages("mvtnorm")
install.packages("pracma")
library(mvtnorm)
library(pracma)

M=350
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
var=array(0, c(2, 2, M))
g = array(0,c(2,2,M))
S_true = matrix(c(1,0.6,0.6,1),ncol=2,nrow=2)
I = diag(1,ncol=k,nrow=k)
S_12=numeric(M)
R_m=numeric(M)

#Funktion fuer DGP
DGP_funktion = function(n){
  r=matrix(0, ncol=k, nrow=n)
  sigma.sim = array(0, c(k, k, n))
  eps = matrix(0,ncol=n,nrow=k)
  Q = array(0, c(k, k, n))
  R = array(0, c(k, k, n))
  D = array(0, c(k, k, n)); Dinv = array(0, c(k, k, n))
  D[,,1] = matrix(c(sqrt(gamma.sim[1]/(1-alpha.sim[1]-beta.sim[1])),0,0,sqrt(gamma.sim[2]/(1-alpha.sim[2]-beta.sim[2]))),ncol=2,nrow=2) #ge?ndert
  Dinv[,,1] = solve(D[,,1])                
  Q[,,1] = S_true  
  #Datengenerierung
  
  sigma.sim[,,1] = matrix(c(gamma.sim[1]/(1-alpha.sim[1]-beta.sim[1]),0,0,gamma.sim[2]/(1-alpha.sim[2]-beta.sim[2])),ncol=2,nrow=2)
  for (t in 2:n){
    sigma.sim[,,t] = diag(c(gamma.sim[1],gamma.sim[2])) + diag(c(alpha.sim[1],alpha.sim[2]))*r[t-1,]%*%t(r[t-1,]) + diag(c(beta.sim[1],beta.sim[2]))*sigma.sim[,,t-1] 
    D[,,t] = diag(sqrt(c(sigma.sim[1,1,t], sigma.sim[2,2,t])))
    Dinv[,,t] = solve(D[,,t]) 
    I = diag(1,ncol=k,nrow=k)
    Q[,,t] = S_true*(1-A-B) + A*((Q[,,t-1]*I)^0.5%*%Dinv[,,t-1]%*%r[t-1,]%*%t(Dinv[,,t-1]%*%r[t-1,])%*%(Q[,,t-1]*I)^0.5) + B*Q[,,t-1]
    Qast <- diag(c((diag(Q[, , t])^-0.5)),2,2)
    R[, , t] <- Qast %*% Q[, , t] %*% Qast
    r[t,] = rmvnorm(n=1, mean=rep(0,k), sigma=D[,,t]%*%R[,,t]%*%D[,,t]) #evtl NUR auf diagolen abs nehmen?
    eps[,t]=Dinv[,,t-1]%*%r[t-1,]
     }
  
  
  
  result=list(r,R,Q,t(eps)) #1:r#2:R#3:Q#4:eps
  return(result)
}

#Monte Carlo Durchläufe
for (i in 1:M){
DGP=DGP_funktion(n)
R_true=array(c(unlist(DGP[2])),c(k,k,n))
R_true=R_true[1,2,]
r=matrix(c(unlist(DGP[1])),ncol=k,nrow=n)
#k Univariate schaetzungen
for (j in 1:k){
  r1=r[,j]
  sigma2 = vector(length=n)
  LL = function(theta) 
  {
    sigma2[1] = var(r1)                                        ###Startwert
    gamma = abs(theta[1])
    alpha = abs(theta[2])
    beta = abs(theta[3])
    for(u in 2:n)
    {
      
      sigma2[u] = gamma + alpha*r1[u-1]^2 + beta*sigma2[u-1]
    }
    
    L =-(1/2)*sum(log(abs(sigma2)))-sum((1/(2*abs(sigma2)))*(r1^2))
    return(-L)

  }
  schaetzung[i,,j] = optim(c(gamma.sim[j],alpha.sim[j],beta.sim[j]),LL, method="L-BFGS-B",lower=c(1e-6,1e-6,1e-6),upper=c(1-1e-6,1-1e-6,1-1e-6))$par
}

#Matrix D
D = array(0, c(k, k, n))
D[,,1] = matrix(c(sqrt(var(r[,1])),0,0,sqrt(var(r[,2]))),ncol=2,nrow=2) 
sigma=array(0, c(k, k, n))
sigma[,,1] = matrix(c(var(r[,1]),0,0,var(r[,2])),ncol=2,nrow=2)
for (t in 2:n){
  sigma[,,t] = diag(c(schaetzung[i,1,1],schaetzung[i,1,2])) + diag(c(schaetzung[i,2,1],schaetzung[i,2,2]))*r[t-1,]%*%t(r[t-1,]) + diag(c(schaetzung[i,3,1],schaetzung[i,3,2]))*sigma[,,t-1] 
  D[,,t] = diag(sqrt(c(sigma[1,1,t], sigma[2,2,t])))
}
D1=matrix(c(D[1,1,]),ncol=1,nrow=n) 
D2=matrix(c(D[2,2,]),ncol=1,nrow=n) 
u1 = as.matrix(r[,1]/D1)
u2 = as.matrix(r[,2]/D2)
eps = cbind(u1,u2)
#cDCC Schaetzung
T <- nrow(eps)
N <- ncol(eps)
Q <- array(NA, c(N, N, T))
R <- array(NA, c(N, N, T))
dcc.11.logLike=
  function (theta)
  {
    logLikeC <- function(Rt, epst) {
      
      return(as.numeric(-1/2 * (log(det(Rt)) + t(epst) %*% solve(Rt) %*% epst)))
    }
    
    
    alpha <- theta[1]
    beta <- theta[2]
                             
    likelihoods = vector(length = n)
    if (alpha+beta<1) {
      
      # Schleife fuer Q*
      Qdiag = array(0,c(N,N,T))
      Qdiag[,,1] = matrix(c(1,0,0,1),ncol=2,nrow=2)       #Initialwert
      for (b in 1:N){
        for (v in 2:T){
          Qdiag[b,b,v] = (1-alpha-beta)*1 + alpha*(eps[v-1,b]^2)*Qdiag[b,b,v-1] + beta*Qdiag[b,b,v-1]
        }
      }
      
      # Berechnung von S
      S=array(0, c(k, k, 1))
      for (g in 1:T){
        S[,,]=S[,,] + Qdiag[,,g]^0.5%*%eps[g,]%*%t(eps[g,])%*%Qdiag[,,g]^0.5
      }
      S=S/T
      S=matrix(S,ncol=2,nrow=2)
      
      Q[,,1]=S
      for (t in 2:T) {
        Q[, , t] <- S*(1 - alpha - beta) + alpha*((Q[,,t-1]*I)^0.5 %*% eps[t-1,] %*% t(eps[t-1,])%*%(Q[,,t-1]*I)^0.5) + beta * Q[,,t-1]
        Qast <- diag(c((diag(Q[, , t])^-0.5)),2,2) 
        R[, , t] <- Qast %*% Q[, , t] %*% Qast 
        tLike <- logLikeC(R[, , t], eps[t, ]) 
        likelihoods[t] <- tLike
      }
      return(-sum(likelihoods))
    }
    else {
      return(1e06)
    }
  }

dcc=optim(c((A-0.02),(B-0.02)),method = "L-BFGS-B",dcc.11.logLike,lower=c(1e-6,1e-6),upper=c(1-1e-6,1-1e-6),hessian=F)
dcc.schaetz[i,] = dcc$par


# Schleife fuer Q*
Qdiag = array(0,c(k,k,n))
Qdiag[,,1] = matrix(c(1,0,0,1),ncol=2,nrow=2)       #Initialwert
for (b in 1:k){
  for (v in 2:n){
    Qdiag[b,b,v] = (1-dcc.schaetz[i,1]-dcc.schaetz[i,2])*1 + dcc.schaetz[i,1]*(eps[v-1,b]^2)*Qdiag[b,b,v-1] + dcc.schaetz[i,2]*Qdiag[b,b,v-1]
  }
}

#Berechnung von S
S=array(0, c(k, k, 1))
for (g in 1:n){
  S[,,]=S[,,] + (Qdiag[,,g]^0.5)%*%eps[g,]%*%t(eps[g,])%*%(Qdiag[,,g]^0.5)
}
S=matrix(S,ncol=2,nrow=2)/n
S_12[i]=S[1,2]
  
#Berechnung von R
Q[,,1]=S
Qast <- diag(c((diag(Q[, , 1])^-0.5)),2,2)
R[,,1]=Qast %*% Q[, , 1] %*% Qast 
for (t in 2:T) {
  Q[, , t] <- S*(1 - dcc.schaetz[i,1] - dcc.schaetz[i,2]) + dcc.schaetz[i,1]*((Q[,,t-1]*I)^0.5 %*% eps[t-1,] %*% t(eps[t-1,])%*%(Q[,,t-1]*I)^0.5) + dcc.schaetz[i,2] * Q[,,t-1]
  Qast <- diag(c((diag(Q[, , t])^-0.5)),2,2) 
  R[, , t] <- Qast %*% Q[, , t] %*% Qast 

}
R_12=R[1,2,]
#mean von R ueber alle t
R_m[i]=mean(abs(R_12-R_true))
}


delete=numeric(length(dcc.schaetz[,2]))
#Optimierer nicht angelaufen oder nicht sinnvoll = entfernen
for (p in 1:length(dcc.schaetz[,2])){
  if (dcc.schaetz[p,2]<0.1|dcc.schaetz[p,2]==B|dcc.schaetz[p,1]==A|dcc.schaetz[p,1]<0.00001){
    delete[p]=p
  }
}
delete=delete[ delete != 0 ]
if (length(delete)>0){
dcc.schaetz=dcc.schaetz[-delete,]
schaetzung=schaetzung[-delete,,]
R_m=R_m[-delete]
}

#Standardabweichungen der cDCC Parameter
#"Pracma" Package notwendig.
sd.dcc=matrix(0,ncol=2,nrow=length(dcc.schaetz[,2]))
likelihood_ges=function(theta){
  
  alpha.dcc <- abs(theta[1])
  beta.dcc <- abs(theta[2])
  
  gamma1 = abs(theta[3])
  alpha1 = abs(theta[4])
  beta1 = abs(theta[5])
  gamma2 = abs(theta[6])
  alpha2 = abs(theta[7])
  beta2 = abs(theta[8])
  
  D = array(0, c(k, k, n))
  D[,,1] = matrix(c(sqrt(var(r[,1])),0,0,sqrt(var(r[,2]))),ncol=2,nrow=2) 
  sigma=array(0, c(k, k, n))
  sigma[,,1] = matrix(c(var(r[,1]),0,0,var(r[,2])),ncol=2,nrow=2)
  for (t in 2:n){
    sigma[,,t] = diag(c(gamma1,gamma2)) + diag(c(alpha1,alpha2))*r[t-1,]%*%t(r[t-1,]) + diag(c(beta1,beta2))*sigma[,,t-1] 
    D[,,t] = diag(sqrt(c(sigma[1,1,t], sigma[2,2,t])))
  }
  D1=matrix(c(D[1,1,]),ncol=1,nrow=n) 
  D2=matrix(c(D[2,2,]),ncol=1,nrow=n) 
  u1 = as.matrix(r[,1]/D1)
  u2 = as.matrix(r[,2]/D2)
  eps = cbind(u1,u2)
  
  logLike <- function(Rt, epst, Dt, rt) {
    
    return(-1/2 * (n*log(2*pi) + 2*log(det(Dt)) + t(rt)%*%solve(Dt)%*%solve(Dt)%*%rt + log(det(Rt)) + t(epst) %*% solve(Rt) %*% epst - t(epst)%*% epst ))
  }
  
  
  T <- nrow(eps)
  N <- ncol(eps)
  Q <- array(NA, c(N, N, T))
  R <- array(NA, c(N, N, T))
  likelihoods = vector(length = n)
  if (alpha.dcc+beta.dcc<1) {
    
    # Schleife fuer Q*
    Qdiag = array(0,c(N,N,T))
    Qdiag[,,1] = matrix(c(1,0,0,1),ncol=2,nrow=2)       #Initialwert
    for (b in 1:N){
      for (v in 2:T){
        Qdiag[b,b,v] = (1-alpha.dcc-beta.dcc)*1 + alpha.dcc*(eps[v-1,b]^2)*Qdiag[b,b,v-1] + beta.dcc*Qdiag[b,b,v-1]
      }
    }
    
    # Berechnung von S
    S=array(0, c(k, k, 1))
    for (g in 1:T){
      S[,,]=S[,,] + Qdiag[,,g]^0.5%*%eps[g,]%*%t(eps[g,])%*%Qdiag[,,g]^0.5
    }
    S=S/T
    S=matrix(S,ncol=2,nrow=2)
    Q[,,1] = S
    Qast <- diag(c((diag(Q[, , 1])^-0.5)),2,2) 
    R[,,1] = Qast %*% Q[, , t] %*% Qast
    for (t in 2:T) {
      Q[, , t] <- S*(1 - alpha.dcc - beta.dcc) + alpha.dcc*((Q[,,t-1]*I)^0.5 %*% eps[t-1,] %*% t(eps[t-1,])%*%(Q[,,t-1]*I)^0.5) + beta.dcc * Q[,,t-1]
      Qast <- diag(c((diag(Q[, , t])^-0.5)),2,2) 
      R[, , t] <- Qast %*% Q[, , t] %*% Qast 
      tLike <- logLike(R[, , t], eps[t, ], D[,,t], r[t,]) 
      likelihoods[t] <- tLike
    }
    return(-sum(likelihoods))
  }
  else {
    return(1e06)
  }
} 
for (i in 1:length(dcc.schaetz[,2])){
  print(i)
  hesse=hessian(f=likelihood_ges,x0=c(dcc.schaetz[i,1],dcc.schaetz[i,2],schaetzung[i,1,1],schaetzung[i,2,1],schaetzung[i,3,1],schaetzung[i,1,2],schaetzung[i,2,2],schaetzung[i,3,2]))
  cov=solve(hesse)
  sd.dcc[i,]=sqrt(c(cov[1,1],cov[2,2]))                 
}

sd.dcc1=sd.dcc[,1]
sd.dcc2=sd.dcc[,2]
#NA entfernen
sd.dcc1=sd.dcc1[!is.na(sd.dcc1)]
sd.dcc2=sd.dcc2[!is.na(sd.dcc2)]




###Output aus M=300 MC Durchläufen. Dadurch einheitlich vergleichbar.
#mean SE
mean(sd.dcc1[1:300]);mean(sd.dcc2[1:300])
#mean Koef
mean(dcc.schaetz[1:300,1])
mean(dcc.schaetz[1:300,2])
#Bias Koef
mean((dcc.schaetz[1:300,1]-A)/A); mean((dcc.schaetz[1:300,2]-B)/B)
#Bias Korrelation
R_MAE=mean(R_m[1:300])
R_MAE

