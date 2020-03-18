###Vergleich der Mehrschrittprognose Methode 1 vs Methode 2.
install.packages("mvtnorm")
library(mvtnorm)

#leere Matrizen und Vektoren
M=10000
n=1000
gamma.sim=c(0.01,0.01); k=length(gamma.sim)
schaetzung=array(0, c(M,3,k))               
D=array(0, c(k, k, n))
dcc.schaetz=matrix(0,ncol=2,nrow=M)
alpha.sim=c(0.05,0.05)
beta.sim=c(0.9,0.9)
A = 0.01
B = 0.98
sd=matrix(0,ncol=2,nrow=M)
abl1 = matrix(0,ncol=1,nrow=M);abl2 = matrix(0,ncol=1,nrow=M)
var=array(0, c(2, 2, M))
g = array(0,c(2,2,M))
S_true = matrix(c(1,0.8,0.8,1),ncol=2,nrow=2)

#Funktion fuer DGP
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
  
  result=list(r,R,Q,t(eps),sigma.sim) #1:r#2:R#3:Q#4:eps#5:Varianz
  return(result)
}

#Funktion fuer forecast Varianz
fcst_h_function=function(K){
  sigma_fcst=array(0,c(k,k,K))

  sigma_fcst[,,1]= diag(gamma.sim) + diag(alpha.sim)*(r[n-K,]%*%t(r[n-K,])) + diag(beta.sim)*sigma[,,n-K]
  

  for (k in 2:K){
    sum=0
    for (j in 0:(k-2)){
      sum = sum + diag(gamma.sim)*(diag(alpha.sim)+diag(beta.sim))^j  
    }
    sigma_fcst[,,k] = sum + ((diag(alpha.sim)+diag(beta.sim)))^(k-1)*sigma_fcst[,,1]
  }
  return(sigma_fcst)
}

#Funktion fuer foreacst Methode 1 
fcst1_R_function=function(K){
  Q_fcst=array(0,c(2,2,K))
  R_fcst1=array(0,c(2,2,K))
  
  Q_fcst[,,1] = (1-A-B)*S_true + A*(eps[n-K,]%*%t(eps[n-K,])) + B*Q[,,n-K] 
  Qast = diag(c((diag(Q_fcst[,,1])^-0.5)),2,2) 
  R_fcst1[, , 1] <- Qast%*%Q_fcst[,,1]%*%Qast 
  
  for (k in 2:K){
    sum=0
    for (j in 0:(k-2)){
      sum = sum + (1-A-B)*S_true*(A+B)^j  
    }
    Q_fcst[,,k] = sum + ((A+B)^(k-1))*Q_fcst[,,1]
    Qast = diag(c((diag(Q_fcst[,,k])^-0.5)),2,2) 
    R_fcst1[, , k] <- Qast %*% Q_fcst[, , k] %*% Qast 
  }
  return(R_fcst1)
}
#Funktion fuer foreacst Methode 2 
fcst2_R_function=function(K){
  R_fcst2=array(0,c(k,k,K))
  Q_fcst = (1-A-B)*S_true + A*(eps[n-K,]%*%t(eps[n-K,])) + B*Q[,,n-K] 
  Qast = diag(c((diag(Q_fcst)^-0.5)),2,2) 
  R_fcst2[, , 1] <- Qast%*%Q_fcst%*%Qast 
  #In Simulation kein Unterschied, da wahres S genutzt wird.
  #Sast=diag(c((diag(S_true)^-0.5)),2,2) 
  #R_quer= Sast%*%S_true%*%Sast      
  R_quer=S_true                    
  for (k in 2:K){
    sum=0
    for (j in 0:(k-2)){
      sum = sum + (1-A-B)*R_quer*(A+B)^j  
    }
    R_fcst2[,,k]= sum + ((A+B)^(k-1))*R_fcst2[,,1]   
  }
  return(R_fcst2)
}


#Monte Carlo Forecast Vergleich Methode 1 und Methode 2
K=250
fcst1_R_collect=matrix(0,ncol=M,nrow=K)
fcst2_R_collect=matrix(0,ncol=M,nrow=K)
R_collect=matrix(0,ncol=M,nrow=n)
for (i in 1:M){
  print(i)
  #DGP
  DGP=DGP_funktion(n)
  sigma=array(c(unlist(DGP[5])),c(k,k,n))
  r=matrix(c(unlist(DGP[1])),ncol=k,nrow=n)
  eps=matrix(c(unlist(DGP[4])),ncol=k,nrow=n) 
  Q=array(c(unlist(DGP[3])),c(k,k,n))
  R=array(c(unlist(DGP[2])),c(k,k,n))
  R_collect[,i]=R[1,2,]
  #Methode 1
  fcst1_R=fcst1_R_function(K=250)
  fcst1_R_collect[,i]=fcst1_R[1,2,]       
  #Methode 2
  fcst2_R=fcst2_R_function(K=250)
  fcst2_R_collect[,i]=fcst2_R[1,2,]
  #Varianz
  fcst_h=fcst_h_function(K=250)
}
R_mean=rowMeans(R_collect)
#Bias Methode 1
fcst1_R_mean=rowMeans(fcst1_R_collect)
bias1=fcst1_R_mean-R_mean[(n-K+1):n]
#Bias Methode 2
fcst2_R_mean=rowMeans(fcst2_R_collect)
bias2=fcst2_R_mean-R_mean[(n-K+1):n]


#Plot 
par(mar=c(5,7,5,5))
plot(bias1,type="l", xlab = "Horizont k", ylab="Bias",main="Korrelation = 0.8",lwd=2,cex.lab=1.8,cex.main=1.9,cex.axis=1.1)
lines(bias2,type="l",col="blue",lwd=2)

