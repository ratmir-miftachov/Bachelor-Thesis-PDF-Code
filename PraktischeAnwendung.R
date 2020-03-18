###packages
install.packages("tseries")
install.packages("rgenoud")
install.packages("forecast")
install.packages("fGarch")
install.packages("zoo")
install.packages("ggplot2")
install.packages("xts")
install.packages("pracma")
library("tseries")
library("rgenoud")
library("forecast")
library("fGarch")
library("zoo")
library("ggplot2")
library("xts")
library("pracma")


#Datenanalyse

#1.Dateneinlesen
#Dax (GER)
#daten_dax = read.zoo("C:/Users/mifra/Desktop/THESIS/ThesisCode/Daten/DAX.csv", sep=",", header=T)
daten_dax = read.zoo(file.choose(), sep=",", header=T)
#SP500 (USA)
#daten_SP=read.zoo("C:/Users/mifra/Desktop/THESIS/ThesisCode/Daten/SP.csv", sep=",", header=T)
daten_SP=read.zoo(file.choose(), sep=",", header=T)
#IBM 
#daten_IBM=read.zoo("C:/Users/mifra/Desktop/THESIS/ThesisCode/Daten/IBM.csv", sep=",", header=T)
daten_IBM=read.zoo(file.choose(), sep=",", header=T)
#IBEX(Spanien)
#daten_IBEX=read.zoo("C:/Users/mifra/Desktop/THESIS/ThesisCode/Daten/IBEX.csv", sep=",", header=T)
daten_IBEX=read.zoo(file.choose(), sep=",", header=T)
#NIKKEI(Japan)
#daten_nik=read.zoo("C:/Users/mifra/Desktop/THESIS/ThesisCode/Daten/NIKKEI.csv", sep=",", header=T)
daten_nik=read.zoo(file.choose(), sep=",", header=T)


###Zusammenfuegen zu einer Matrix
data=cbind(daten_dax$Adj.Close,daten_SP$Adj.Close,daten_IBM$Adj.Close,daten_IBEX$Adj.Close,daten_nik$Adj.Close)
data=zoo(na.omit(data))  
time=index(data)
time=time[-1]
kurs=matrix(data,ncol=5)
r=diff(log(kurs))

#Spalten- und Zeilenlaengen
k=ncol(r)
n=length(r[,1])
#Falls Rendite=0, dann ersetze durch Rendite Vortag. 
for (j in 1:k){
  for (i in 1:n){
    if (r[i,j]==0){
      r[i,j]=r[i-1,j]
      
    }
  }
}


#leere Matrizen und Vektoren
schaetzung=matrix(0, ncol=3,nrow=k) 
sd=matrix(0, ncol=3,nrow=k) 
I = diag(1,ncol=k,nrow=k)
dcc.schaetz=matrix(0,ncol=k,nrow=1)
se_uni=array(0,c(3,3,5))
#p-Wert
t1=numeric(k)
t2=numeric(k)
t3=numeric(k)
p.value1=numeric(k)
p.value2=numeric(k)
p.value3=numeric(k)



# 1.Schritt: k Univariate GARCH Schaetzungen
  for (j in 1:k){
    r1=r[,j]
    sigma2 = vector(length=n)
    LL = function(theta) 
    {
      sigma2[1] = var(r1)                                       
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
    garch = nlminb(c(0.00001,0.065,0.93),LL,lower=c(1e-6,1e-6,1e-6),upper=c(1-1e-6,1-1e-6,1-1e-6))
    schaetzung[j,]=garch$par
    Hesse=optimHess(c(schaetzung[j,]),LL)
    cov=solve(Hesse)
    sd[j,]=sqrt(diag(cov))
    t1[j]=schaetzung[j,1]/sd[j,1]
    t2[j]=schaetzung[j,2]/sd[j,2]
    t3[j]=schaetzung[j,3]/sd[j,3]
    p.value1[j] = 2*pt(abs(t1[j]), df=n-(2+1),lower=FALSE)
    p.value2[j] = 2*pt(abs(t2[j]), df=n-(2+1),lower=FALSE)
    p.value3[j] = 2*pt(abs(t3[j]), df=n-(2+1),lower=FALSE)
        }

###output GARCH
#GARCH Schaetzer
schaetzung
#garch Standardabweichungen
sd
#garch p-Werte
p.value1
p.value2
p.value3


# 2.Schritt: Renditen standardisieren
D = array(0, c(k, k, n))
D[,,1] = diag(c(sqrt(var(r[,1])),sqrt(var(r[,2])),sqrt(var(r[,3])),sqrt(var(r[,4])),sqrt(var(r[,5]))))  
sigma=array(0, c(k, k, n))
sigma[,,1] =diag(c(var(r[,1]),var(r[,2]),var(r[,3]),var(r[,4]),var(r[,5])))
for (t in 2:n){
  sigma[,,t] = diag(c(schaetzung[,1])) + diag(c(schaetzung[,2]))*r[t-1,]%*%t(r[t-1,]) + diag(c(schaetzung[,3]))*sigma[,,t-1] 
  D[,,t] = sqrt(sigma[,,t])
}
D1=matrix(c(D[1,1,]),ncol=1,nrow=n) 
D2=matrix(c(D[2,2,]),ncol=1,nrow=n) 
D3=matrix(c(D[3,3,]),ncol=1,nrow=n)
D4=matrix(c(D[4,4,]),ncol=1,nrow=n) 
D5=matrix(c(D[5,5,]),ncol=1,nrow=n) 
u1 = as.matrix(r[,1]/D1)
u2 = as.matrix(r[,2]/D2)
u3 = as.matrix(r[,3]/D2)
u4 = as.matrix(r[,4]/D2)
u5 = as.matrix(r[,5]/D2)
eps = cbind(u1,u2,u3,u4,u5)

#3.Schtitt: cDCC Schaetzung
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
        
        #4.1 Schleife fuer Q*
        Qdiag = array(0,c(N,N,T))
        Qdiag[,,1] = diag(1,ncol=5,nrow=5)       #Initialwert
        for (b in 1:N){
          for (v in 2:T){
            Qdiag[b,b,v] = (1-alpha-beta)*1 + alpha*(eps[v-1,b]^2)*Qdiag[b,b,v-1] + beta*Qdiag[b,b,v-1]
          }
        }
        
        #4.2 Berechnung von S
        S=array(0, c(k, k, 1))
        for (g in 1:T){
          S[,,]=S[,,] + Qdiag[,,g]^0.5%*%eps[g,]%*%t(eps[g,])%*%Qdiag[,,g]^0.5
        }
        S=S/T
        S=matrix(S,ncol=5,nrow=5)
        
        Q[,,1]=S
        for (t in 2:T) {
          Q[, , t] <- S*(1 - alpha - beta) + alpha*((Q[,,t-1]*I)^0.5 %*% eps[t-1,] %*% t(eps[t-1,])%*%(Q[,,t-1]*I)^0.5) + beta * Q[,,t-1]
          Qast <- diag(c((diag(Q[, , t])^-0.5)),5,5) 
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
dcc=optim(c(0.05,0.9),dcc.11.logLike,method="L-BFGS-B",lower=c(1e-6,1e-6),upper=c(1-1e-6,1-1e-6))
dcc.schaetz = dcc$par
#cDCC Schaetzer
dcc.schaetz

###Berechnung S fuer Korrelation
# Schleife fuer Q*
Qdiag = array(0,c(N,N,T))
Qdiag[,,1] = diag(1,ncol=5,nrow=5)       
for (b in 1:N){
  for (v in 2:T){
    Qdiag[b,b,v] = (1-dcc.schaetz[1]-dcc.schaetz[2])*1 + dcc.schaetz[1]*(eps[v-1,b]^2)*Qdiag[b,b,v-1] + dcc.schaetz[2]*Qdiag[b,b,v-1]
  }
}
# Berechnung von S
S=array(0, c(k, k, 1))
for (g in 1:T){
  S[,,]=S[,,] + (Qdiag[,,g]^0.5)%*%eps[g,]%*%t(eps[g,])%*%(Qdiag[,,g]^0.5)
}
S=S/T
S=matrix(S,ncol=5,nrow=5)
#Berechnung der Korrelation
Q[,,1]=S
Qast <- diag(c((diag(Q[, , 1])^-0.5)),5,5)
R[,,1]=Qast %*% Q[, , 1] %*% Qast 
for (t in 2:T) {
    Q[, , t] <- S*(1 - dcc.schaetz[1] - dcc.schaetz[2]) + dcc.schaetz[1]*((Q[,,t-1]*I)^0.5 %*% eps[t-1,] %*% t(eps[t-1,])%*%(Q[,,t-1]*I)^0.5) + dcc.schaetz[2] * Q[,,t-1]
    Qast <- diag(c((diag(Q[, , t])^-0.5)),5,5) 
    R[, , t] <- Qast %*% Q[, , t] %*% Qast 
    
  }


#"Pracma" Package notwendig.
#vollständige ML Funktion
likelihood_ges=function(theta){
  #1.Univariat
  
  
  alpha.dcc <- abs(theta[1])
  beta.dcc <- abs(theta[2])
  
  gamma1 = abs(theta[3])
  gamma2 = abs(theta[4])
  gamma3 = abs(theta[5])
  gamma4 = abs(theta[6])
  gamma5 = abs(theta[7])
  alpha1 = abs(theta[8])
  alpha2 = abs(theta[9])
  alpha3 = abs(theta[10])
  alpha4 = abs(theta[11])
  alpha5 = abs(theta[12])
  beta1 =  abs(theta[13])
  beta2 =  abs(theta[14])
  beta3 =  abs(theta[15])
  beta4 =  abs(theta[16])
  beta5 =  abs(theta[17])

  
  #2.eps
  D = array(0, c(k, k, n))
  D[,,1] = diag( c(sqrt(var(r[,1])), sqrt(var(r[,2])), sqrt(var(r[,3])), sqrt(var(r[,4])), sqrt(var(r[,5]))),ncol=5,nrow=5 ) 
  sigma=array(0, c(k, k, n))
  sigma[,,1] = diag(c(var(r[,1]),var(r[,2]),var(r[,3]),var(r[,4]),var(r[,5])),ncol=5,nrow=5)
  for (t in 2:n){
    sigma[,,t] = diag(c(gamma1,gamma2,gamma3,gamma4,gamma5)) + diag(c(alpha1,alpha2,alpha3,alpha4,alpha5))*r[t-1,]%*%t(r[t-1,]) + diag(c(beta1,beta2,beta3,beta4,beta5))*sigma[,,t-1] 
    D[,,t] = sqrt(sigma[,,t])
  }
  D1=matrix(c(D[1,1,]),ncol=1,nrow=n) 
  D2=matrix(c(D[2,2,]),ncol=1,nrow=n) 
  D3=matrix(c(D[3,3,]),ncol=1,nrow=n)
  D4=matrix(c(D[4,4,]),ncol=1,nrow=n) 
  D5=matrix(c(D[5,5,]),ncol=1,nrow=n) 
  u1 = as.matrix(r[,1]/D1)
  u2 = as.matrix(r[,2]/D2)
  u3 = as.matrix(r[,3]/D2)
  u4 = as.matrix(r[,4]/D2)
  u5 = as.matrix(r[,5]/D2)
  eps = cbind(u1,u2,u3,u4,u5)
  #3.DCC
  logLike <- function(Rt, epst, Dt, rt) {
    
    return(-1/2 * (n*log(2*pi) + 2*log(det(Dt)) + t(rt)%*%solve(Dt)%*%solve(Dt)%*%rt + log(det(Rt)) + t(epst) %*% solve(Rt) %*% epst - t(epst)%*% epst ))
  }
  
  T <- nrow(eps)
  N <- ncol(eps)
  Q <- array(NA, c(N, N, T))
  R <- array(NA, c(N, N, T))
  Q[, , 1] <- S 
  likelihoods = vector(length = n)
  if (alpha.dcc+beta.dcc<1) {
    
    #4.1 Schleife fuer Q*
    Qdiag = array(0,c(N,N,T))
    Qdiag[,,1] = diag(1,ncol=5,nrow=5)      #Initialwert
    for (b in 1:N){
      for (v in 2:T){
        Qdiag[b,b,v] = (1-alpha.dcc-beta.dcc)*1 + alpha.dcc*(eps[v-1,b]^2)*Qdiag[b,b,v-1] + beta.dcc*Qdiag[b,b,v-1]
      }
    }
    
    #4.2 Berechnung von S
    S=array(0, c(k, k, 1))
    for (g in 1:T){
      S[,,]=S[,,] + Qdiag[,,g]^0.5%*%eps[g,]%*%t(eps[g,])%*%Qdiag[,,g]^0.5
    }
    S=S/T
    S=matrix(S,ncol=5,nrow=5)
    Q[,,1] = S
    Qast <- diag(c((diag(Q[, , 1])^-0.5)),5,5) 
    R[,,1] = Qast %*% Q[, , t] %*% Qast
    for (t in 2:T) {
      Q[, , t] <- S*(1 - alpha.dcc - beta.dcc) + alpha.dcc*((Q[,,t-1]*I)^0.5 %*% eps[t-1,] %*% t(eps[t-1,])%*%(Q[,,t-1]*I)^0.5) + beta.dcc * Q[,,t-1]
      Qast <- diag(c((diag(Q[, , t])^-0.5)),5,5) 
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
#Standardabweichungen der cDCC Parameter bestimmen.
hesse=hessian(f=likelihood_ges,x0=c(dcc.schaetz[1],dcc.schaetz[2],schaetzung[1,1],schaetzung[2,1],schaetzung[3,1],schaetzung[4,1],schaetzung[5,1],schaetzung[1,2],schaetzung[2,2],schaetzung[3,2],schaetzung[4,2],schaetzung[5,2],schaetzung[1,3],schaetzung[2,3],schaetzung[3,3],schaetzung[4,3],schaetzung[5,3]))
cov=solve(hesse)
sd.dcc=sqrt(c(cov[1,1],cov[2,2]))                 

t1.dcc=dcc.schaetz[1]/sd.dcc[1]
t2.dcc=dcc.schaetz[2]/sd.dcc[2]
p.value1.dcc = 2*pt(abs(t1.dcc), df=n-(2+1),lower=FALSE)
p.value2.dcc = 2*pt(abs(t2.dcc), df=n-(2+1),lower=FALSE)

#Standardabweichung cDCC Parameter
sd.dcc
#p-Werte cdcc Parameter
p.value1.dcc
p.value2.dcc

###Deskriptive Statistik
summary(r)
kurtosis(r[,1],method="excess")
kurtosis(r[,2],method="excess")
kurtosis(r[,3],method="excess")
kurtosis(r[,4],method="excess")
kurtosis(r[,5],method="excess")


#PLOTS

#Renditen Plots
par(mar=c(5,7,5,5))
plot(r[,1]~time,type="l",xlab="Zeit",ylab="Log Renditen",main="DAX 30",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(r[,2]~time,type="l",xlab="Zeit",ylab="Log Renditen",main="S&P 500",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(r[,3]~time,type="l",xlab="Zeit",ylab="Log Renditen",main="IBM",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(r[,4]~time,type="l",xlab="Zeit",ylab="Log Renditen",main="IBEX 35",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(r[,5]~time,type="l",xlab="Zeit",ylab="Log Renditen",main="Nikkei 225",cex.lab=2.7,cex.main=3,cex.axis=2)
#Volatilitaeten Plots
plot(D[1,1,]~time,type="l",xlab="Zeit",ylab="VolatilitÃ¤ten",main="DAX 30",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(D[2,2,]~time,type="l",xlab="Zeit",ylab="VolatilitÃ¤ten",main="S&P 500",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(D[3,3,]~time,type="l",xlab="Zeit",ylab="VolatilitÃ¤ten",main="IBM",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(D[4,4,]~time,type="l",xlab="Zeit",ylab="VolatilitÃ¤ten",main="IBEX 35",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(D[5,5,]~time,type="l",xlab="Zeit",ylab="VolatilitÃ¤ten",main="Nikkei 225",cex.lab=2.7,cex.main=3,cex.axis=2)
#Korrelationen Plots
plot(R[1,2,]~time,type="l",xlab="Zeit",ylab="Korrelationen",main="DAX vs S&P 500",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(R[1,3,]~time,type="l",xlab="Zeit",ylab="Korrelationen",main="DAX vs IBM",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(R[1,4,]~time,type="l",xlab="Zeit",ylab="Korrelationen",main="DAX vs IBEX 35",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(R[1,5,]~time,type="l",xlab="Zeit",ylab="Korrelationen",main="DAX vs Nikkei 225",cex.lab=2.7,cex.main=3,cex.axis=2)

plot(R[2,3,]~time,type="l",xlab="Zeit",ylab="Korrelationen",main="S&P 500 vs IBM",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(R[2,4,]~time,type="l",xlab="Zeit",ylab="Korrelationen",main="S&P 500 vs IBEX",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(R[2,5,]~time,type="l",xlab="Zeit",ylab="Korrelationen",main="S&P 500 vs Nikkei 225",cex.lab=2.7,cex.main=3,cex.axis=2)

plot(R[3,4,]~time,type="l",xlab="Zeit",ylab="Korrelationen",main="IBM vs IBEX",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(R[3,5,]~time,type="l",xlab="Zeit",ylab="Korrelationen",main="IBM vs Nikkei 225",cex.lab=2.7,cex.main=3,cex.axis=2)

plot(R[4,5,]~time,type="l",xlab="Zeit",ylab="Korrelationen",main="IBEX vs Nikkei 225",cex.lab=2.7,cex.main=3,cex.axis=2)

#ACF der Renditen und der Renditen^2
par(cex.main=3)
plot(acf(r[,1],plot = F),main="DAX",lwd=2,cex.lab=2.7,cex.axis=2)
plot(acf(r[,2],plot = F),main="S&P 500",lwd=2,cex.lab=2.7,cex.axis=2)
plot(acf(r[,3],plot = F),main="IBM",lwd=2,cex.lab=2.7,cex.axis=2)
plot(acf(r[,4],plot = F),main="IBEX",lwd=2,cex.lab=2.7,cex.axis=2)
plot(acf(r[,5],plot = F),main="Nikkei",lwd=2,cex.lab=2.7,cex.axis=2)


plot(acf(r[,1]^2,plot = F),main="DAX",lwd=2,cex.lab=2.7,cex.axis=2)
plot(acf(r[,2]^2,plot = F),main="S&P 500",lwd=2,cex.lab=2.7,cex.axis=2)
plot(acf(r[,3]^2,plot = F),main="IBM",lwd=2,cex.lab=2.7,cex.axis=2)
plot(acf(r[,4]^2,plot = F),main="IBEX",lwd=2,cex.lab=2.7,cex.axis=2)
plot(acf(r[,5]^2,plot = F),main="Nikkei",lwd=2,cex.lab=2.7,cex.axis=2)

################################################################################################
#Mehrschrittprognose mit Methode 2 
S_true = cov(eps)
#Funktion für forecast Varianz
fcst_h_function=function(K){
  sigma_fcst=array(0,c(k,k,K))
  
  sigma_fcst[,,1]= diag(schaetzung[,1]) + diag(schaetzung[,2])*(r[n,]%*%t(r[n,])) + diag(schaetzung[,3])*sigma[,,n]
  
  
  for (k in 2:K){
    sum=0
    for (j in 0:(k-2)){
      sum = sum + diag(schaetzung[,1])*(diag(schaetzung[,2])+diag(schaetzung[,3]))^j  
    }
    sigma_fcst[,,k] = sum + ((diag(schaetzung[,2])+diag(schaetzung[,3])))^(k-1)*sigma_fcst[,,1]
  }
  return(sigma_fcst)
}
#Funktion fuer foreacst Methode 2 
fcst2_R_function=function(K){
  R_fcst2=array(0,c(k,k,K))
  #Rekursion für t+1.
  Q_fcst1 = (1-dcc.schaetz[1]-dcc.schaetz[2])*S_true + dcc.schaetz[1]*(eps[n,]%*%t(eps[n,])) + dcc.schaetz[2]*Q[,,n] 
  Qast = diag(c((diag(Q_fcst1)^-0.5)),5,5) 
  #Reskalierung für Korrelation in t+1
  R_fcst2[,,1]=Qast%*%Q_fcst1%*%Qast
  
  Sast=diag(c((diag(S_true)^-0.5)),5,5)
  R_quer=Sast%*%S_true%*%Sast   #R_quer= cor(eps)

  for (k in 2:K){
    sum=0
    for (j in 0:(k-2)){
      sum = sum + (1-dcc.schaetz[1]-dcc.schaetz[2])*R_quer*(dcc.schaetz[1]+dcc.schaetz[2])^j  
    }
    R_fcst2[,,k]= sum + ((dcc.schaetz[1]+dcc.schaetz[2])^(k-1))*R_fcst2[,,1]   
  }
  return(R_fcst2)
}


#Varianz
fcst_h=sqrt(fcst_h_function(K=364))
#Korrelationen
fcst_R=fcst2_R_function(364)


#n-50 SD
sd_fcst1=D[1,1,(n-49):n]
sd_fcst2=D[2,2,(n-49):n]
sd_fcst3=D[3,3,(n-49):n]
sd_fcst4=D[4,4,(n-49):n]
sd_fcst5=D[5,5,(n-49):n]

sd_fcst11=c(sd_fcst1,fcst_h[1,1,])
sd_fcst22=c(sd_fcst2,fcst_h[2,2,])
sd_fcst33=c(sd_fcst3,fcst_h[3,3,])
sd_fcst44=c(sd_fcst4,fcst_h[4,4,])
sd_fcst55=c(sd_fcst5,fcst_h[5,5,])



#Plots der prog. Standardabweichungen

time_fcst=seq.Date(as.Date("2017-08-23"), as.Date("2018-10-10"),by='days')  
par(mar=c(5,7,5,5))
plot(sd_fcst11~time_fcst,type="l",main="DAX",ylab="Volatilität",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(sd_fcst22~time_fcst,type="l",main="S&P 500",ylab="Volatilität",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(sd_fcst33~time_fcst,type="l",main="IBM",ylab="Volatilität",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(sd_fcst44~time_fcst,type="l",main="IBEX 35",ylab="Volatilität",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(sd_fcst55~time_fcst,type="l",main="Nikkei 225",ylab="Volatilität",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)



##Korrelationen
R_fcst12=R[1,2,(n-49):n]
R_fcst13=R[1,3,(n-49):n]
R_fcst14=R[1,4,(n-49):n]
R_fcst15=R[1,5,(n-49):n]
R_fcst23=R[2,3,(n-49):n]
R_fcst24=R[2,4,(n-49):n]
R_fcst25=R[2,5,(n-49):n]
R_fcst34=R[3,4,(n-49):n]
R_fcst35=R[3,5,(n-49):n]
R_fcst45=R[4,5,(n-49):n]

R_fcst12=c(R_fcst12,fcst_R[1,2,])
R_fcst13=c(R_fcst13,fcst_R[1,3,])
R_fcst14=c(R_fcst14,fcst_R[1,4,])
R_fcst15=c(R_fcst15,fcst_R[1,5,])
R_fcst23=c(R_fcst23,fcst_R[2,3,])
R_fcst24=c(R_fcst24,fcst_R[2,4,])
R_fcst25=c(R_fcst25,fcst_R[2,5,])
R_fcst34=c(R_fcst34,fcst_R[3,4,])
R_fcst35=c(R_fcst35,fcst_R[3,5,])
R_fcst45=c(R_fcst45,fcst_R[4,5,])

#Plots der prog. Korrelationen

plot(R_fcst12~time_fcst,type="l",main="DAX vs S&P 500",ylab="Korrelation",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(R_fcst13~time_fcst,type="l",main="DAX vs IBM",ylab="Korrelation",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(R_fcst14~time_fcst,type="l",main="DAX vs IBEX 35",ylab="Korrelation",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(R_fcst15~time_fcst,type="l",main="DAX vs Nikkei 225",ylab="Korrelation",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)

plot(R_fcst23~time_fcst,type="l",main="S&P 500 vs IBM",ylab="Korrelation",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(R_fcst24~time_fcst,type="l",main="S&P 500 vs IBEX",ylab="Korrelation",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(R_fcst25~time_fcst,type="l",main="S&P 500 vs Nikkei 225",ylab="Korrelation",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)

plot(R_fcst34~time_fcst,type="l",main="IBM vs IBEX",ylab="Korrelation",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)
plot(R_fcst35~time_fcst,type="l",main="IBM vs Nikkei 225",ylab="Korrelation",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)

plot(R_fcst45~time_fcst,type="l",main="IBEX vs Nikkei 225",ylab="Korrelation",xlab="Zeit",cex.lab=2.7,cex.main=3,cex.axis=2)

#################################################################################################

#Aielli VaR forecast (Vorsicht: überschreibt oberen Code teilweise.)
k_hzt=900
T_cut=n-k_hzt
dcc.schaetz=numeric(k_hzt); VaR=numeric(k_hzt);VaR=numeric(k_hzt); sigma_fcst=array(0,c(k,k,k_hzt));D_fcst=array(0,c(k,k,k_hzt));Q_fcst=array(0,c(k,k,k_hzt));R_fcst=array(0,c(k,k,k_hzt));var_fcst=array(0,c(k,k,k_hzt))
k_hzt=k_hzt-1
schaetzung=matrix(0, ncol=3,nrow=k) 
I = diag(1,ncol=k,nrow=k)
w=c(0.2,0.2,0.2,0.2,0.2)

for (c in 0:(k_hzt)){
  T_fcst=T_cut+c
  r_cut=r[1:T_fcst,]                #Daten kuerzen
  
  #Univariate Schaetzungen
  for (j in 1:k){
    r1=r_cut[,j]
    sigma2 = vector(length=T_fcst)
    LL = function(theta) 
    {
      sigma2[1] = var(r1)                                     
      gamma = abs(theta[1])
      alpha = abs(theta[2])
      beta = abs(theta[3])
      for(u in 2:T_fcst)
      {
        
        sigma2[u] = gamma + alpha*r1[u-1]^2 + beta*sigma2[u-1]
      }
      
      L =-(1/2)*sum(log(abs(sigma2)))-sum((1/(2*abs(sigma2)))*(r1^2))
      return(-L)
      
    }
    schaetzung[j,] = optim(c(0.00001,0.065,0.93),LL, method="L-BFGS-B",lower=c(1e-6,1e-6,1e-6),upper=c(1-1e-6,1-1e-6,1-1e-6))$par
  }
  
  D = array(0, c(k, k, T_fcst))
  D[,,1] = diag(sqrt(c(var(r_cut[,1]),var(r_cut[,2]),var(r_cut[,3]),var(r_cut[,4]),var(r_cut[,5])))) 
  sigma=array(0, c(k, k, T_fcst))
  sigma[,,1] = diag(c(var(r_cut[,1]),var(r_cut[,2]),var(r_cut[,3]),var(r_cut[,4]),var(r_cut[,5])))
  for (t in 2:T_fcst){
    sigma[,,t] = diag(c(schaetzung[,1])) + diag(c(schaetzung[,2]))*r_cut[t-1,]%*%t(r_cut[t-1,]) + diag(c(schaetzung[,3]))*sigma[,,t-1] 
    D[,,t] = sqrt(sigma[,,t])
  }
  D1=matrix(c(D[1,1,]),ncol=1,nrow=T_fcst) 
  D2=matrix(c(D[2,2,]),ncol=1,nrow=T_fcst) 
  D3=matrix(c(D[3,3,]),ncol=1,nrow=T_fcst) 
  D4=matrix(c(D[4,4,]),ncol=1,nrow=T_fcst) 
  D5=matrix(c(D[5,5,]),ncol=1,nrow=T_fcst) 
  
  u1 = as.matrix(r_cut[,1]/D1)
  u2 = as.matrix(r_cut[,2]/D2)
  u3 = as.matrix(r_cut[,3]/D1)
  u4 = as.matrix(r_cut[,4]/D1)
  u5 = as.matrix(r_cut[,5]/D1)
  
  eps = cbind(u1,u2,u3,u4,u5)
  
  #cDCC Schaetzung
  dcc.11.logLike=
    function (theta)
    {
      logLikeC <- function(Rt, epst) {
        
        return(as.numeric(-1/2 * (log(det(Rt)) + t(epst) %*% solve(Rt) %*% epst)))
      }
      
      
      alpha <- theta[1]
      beta <- theta[2]
      T <- nrow(eps)
      N <- ncol(eps)
      Q <- array(NA, c(N, N, T))
      R <- array(NA, c(N, N, T))
      likelihoods = vector(length = T_fcst)
      if (alpha+beta<1) {
        
        #4.1 Schleife fÃ¼r Q*
        Qdiag = array(0,c(N,N,T))
        Qdiag[,,1] = diag(1,ncol=5,nrow=5)       #Initialwert
        for (b in 1:N){
          for (v in 2:T){
            Qdiag[b,b,v] = (1-alpha-beta)*1 + alpha*(eps[v-1,b]^2)*Qdiag[b,b,v-1] + beta*Qdiag[b,b,v-1]
          }
        }
        
        #4.2 Berechnung von S
        S=array(0, c(k, k, 1))
        for (g in 1:T){
          S[,,]=S[,,] + Qdiag[,,g]^0.5%*%eps[g,]%*%t(eps[g,])%*%Qdiag[,,g]^0.5
        }
        S=S/T
        S=matrix(S,ncol=5,nrow=5)
        
        Q[,,1]=S
        for (t in 2:T) {
          Q[, , t] <- S*(1 - alpha - beta) + alpha*((Q[,,t-1]*I)^0.5 %*% eps[t-1,] %*% t(eps[t-1,])%*%(Q[,,t-1]*I)^0.5) + beta * Q[,,t-1]
          Qast <- diag(c((diag(Q[, , t])^-0.5)),5,5) 
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
  dcc=optim(c(0.05,0.9),dcc.11.logLike,method="L-BFGS-B",lower=c(1e-6,1e-6),upper=c(1-1e-6,1-1e-6))
  dcc.schaetz = dcc$par
  #fcst Univariat
  sigma_fcst[,,c+1]=diag(c(schaetzung[,1])) + diag(c(schaetzung[,2]))*r_cut[T_fcst,]%*%t(r_cut[T_fcst,]) + diag(c(schaetzung[,3]))*sigma[,,T_fcst]
  D_fcst[,,c+1] = sqrt(sigma_fcst[,,c+1])
  
  #fcst cDCC
  Q <- array(NA, c(k, k, T_fcst))
  R <- array(NA, c(k, k, T_fcst))
  
  
  #Berechnung S mit geschaetzten Koef
  Qdiag = array(0,c(k,k,T_fcst))
  Qdiag[,,1] = diag(1,ncol=5,nrow=5)       
  for (b in 1:k){
    for (v in 2:T_fcst){
      Qdiag[b,b,v] = (1-dcc.schaetz[1]-dcc.schaetz[2])*1 + dcc.schaetz[1]*(eps[v-1,b]^2)*Qdiag[b,b,v-1] + dcc.schaetz[2]*Qdiag[b,b,v-1]
    }
  }
  S=array(0, c(k, k, 1))
  for (g in 1:T_fcst){
    S[,,]=S[,,] + (Qdiag[,,g]^0.5)%*%eps[g,]%*%t(eps[g,])%*%(Qdiag[,,g]^0.5)
  }
  S=S/T_fcst
  S=matrix(S,ncol=5,nrow=5)
  Q[,,1]=S
  Qast <- diag(c((diag(Q[, , 1])^-0.5)),5,5)
  R[,,1]=Qast %*% Q[, , 1] %*% Qast 
  for (t in 2:T_fcst) {
    Q[, , t] <- S*(1 - dcc.schaetz[1] - dcc.schaetz[2]) + dcc.schaetz[1]*((Q[,,t-1]*I)^0.5) %*% eps[t-1,] %*% t(eps[t-1,])%*%((Q[,,t-1]*I)^0.5) + dcc.schaetz[2] * Q[,,t-1]
    Qast <- diag(c((diag(Q[, , t])^-0.5)),5,5) 
    R[, , t] <- Qast %*% Q[, , t] %*% Qast 
    
  }
  
  
  Q_fcst[,,c+1] = S*(1-dcc.schaetz[1]-dcc.schaetz[2]) + dcc.schaetz[1]*(((Q[,,T_fcst]*I)^0.5)%*%eps[T_fcst,]%*%t(eps[T_fcst,])%*%((Q[,,T_fcst]*I)^0.5)) + dcc.schaetz[2]*Q[,,T_fcst]
  Qast <- diag(c((diag(Q_fcst[, , c+1])^-0.5)),5,5) 
  R_fcst[,,c+1] =  Qast%*%Q_fcst[,,c+1]%*%Qast 
  
  var_fcst[,,c+1]=D_fcst[,,c+1]%*%R_fcst[,,c+1]%*%D_fcst[,,c+1]
  rp_var_fcst =t(w)%*%var_fcst[,,c+1]%*%w                
  VaR[c+1]=qnorm(0.05,mean=0,sd=sqrt(rp_var_fcst))
  
}

#UC-Test
r_p=r[(T_cut+1):n,]%*%w
X=sum(r_p<VaR)
alpha_hut=X/length(VaR)
LR=2*(log((alpha_hut^X)*(1-alpha_hut)^(length(VaR)-X)) - log((0.05^X)*(1-0.05)^(length(VaR)-X))) 
abl_rate = LR > qchisq(1-0.05,df=1)     

#VaR plot
par(mar=c(5,7,5,5))
time_cut=time[(T_cut+1):n]
plot(r_p~time_cut,type="l",xlab="Zeit",ylab="Log Renditen",main="Gleichgewichtiges Portfolio",lwd=1.5,cex.lab=1.8,cex.main=1.9,cex.axis=1.1)
lines(VaR~time_cut,col="blue",lwd=1.5)


