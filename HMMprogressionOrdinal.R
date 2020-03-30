##################################################
### Paper:
### A hidden Markov model to address measurement errors in ordinal response scale and non-decreasing process
###
### Authors:
### Lizbeth Naranjo (1), Luz Judith R. Esparza (2), Carlos J. Perez (3).
###
### (1) Departamento de Matemáticas, Facultad de Ciencias, Universidad Nacional Autonoma de Mexico (UNAM), Mexico
### (2) Catedra CONACyT - Universidad Autonoma de Aguascalientes, Mexico
### (3) Departamento de Matematicas, Facultad de Veterinaria, Universidad de Extremadura, Spain
### 
### Journal: Mathematics
### Submitted. Under Revision. 
##################################################

##################################################
### R packages
###
### Load the R libraries
##################################################
library(rjags)
#library(mnormt)
library(MCMCpack) ### MCMC
#library(LearnBayes) ### Bayes
#library(gamlss.dist) ### Inverse Gaussian distribution
##################################################


##################################################
### ADDRESS
###
### Instructions: 
### Change the address where the data and codes are located. 
### setwd("HERE")
##################################################
setwd("~/Documents/Articulos/ArticuloHMMordinal/HMMordinalCodes/")
getwd()
##################################################


##################################################
### SIMULATE DATA
###
### Run the followin lines of code.
##################################################

##################################################
### Truncated Normal I[tra < Z]
rnormleft <- function(tra,mu,sig){
  rp <- pnorm(tra, mean=mu, sd=sig)
  u <- rp+(1-rp)*runif(1)
  q <- qnorm(u, mean=mu, sd=sig)
  if(!is.finite(q)){ q = tra }
  return(q)
}
##################################################

##################################################
### Parameters 

N = 100   # subjects
K = 6   # times
J = 4   # categories

X = array( runif(N*2) ,dim=c(N,2)) 
Z = array( runif(N*K*2) ,dim=c(N,K,2)) 

Beta = c(-1,1)
Gama = c(-1,1)
Kappa = c(-Inf,0,2,4,Inf)
Sigma2 = c(0.5,1)^2

Xi = c(rep(1,N/2),rep(2,N/2))

##################################################
### Simulate Response

par(mfrow=c(1,3))

set.seed(12345)

eta = array(NA,dim=c(N,K))
for(k in 1:K){
  eta[,k] = X[,]%*%Beta + Z[,k,]%*%Gama 
}
Y = W = array(NA,dim=c(N,K))
Ymis = Wmis = array(NA,dim=c(N,K))
for(i in 1:N){
  W[i,1] = rnorm(1,eta[i,1],1)
  Wmis[i,1] = rnorm(1,W[i,1],sqrt(Sigma2[Xi[i]]))
  for(k in 2:K){
    W[i,k] = rnormleft(W[i,k-1],eta[i,k],1)
    Wmis[i,k] = rnorm(1,W[i,k],sqrt(Sigma2[Xi[i]]))
}	}
for(i in 1:N){
  for(k in 1:K){
    for(j in 1:J){
      if(Kappa[j]< Wmis[i,k] & Wmis[i,k]<=Kappa[j+1]){
        Ymis[i,k] = j
      }
      if(Kappa[j]< W[i,k] & W[i,k]<=Kappa[j+1]){
        Y[i,k] = j
} } } }	

plot(W,Wmis)
lines(c(-5,20),c(-5,20),col=2,lwd=2)   
plot(c(1:K),W[1,], main="W", type="l",col=1, ylim=c(min(W,na.rm=T),max(W,na.rm=T)))
for(i in 2:N){	lines(c(1:K),W[i,],col=i)		}
plot(c(1:K),Wmis[1,], main="W mis", type="l",col=1, ylim=c(min(Wmis,na.rm=T),max(Wmis,na.rm=T)))
for(i in 2:N){	lines(c(1:K),Wmis[i,],col=i)	}
  
for(k in 2:K){ print(table(Y[,k-1],Y[,k])) }
for(k in 2:K){ print(table(Ymis[,k-1],Ymis[,k])) }

tabla = tablamis = matrix(0,J,J)
for(i in 1:N){
  for(k in 2:K){ 
    tabla[Y[i,k-1],Y[i,k]] = tabla[Y[i,k-1],Y[i,k]] +1
    tablamis[Ymis[i,k-1],Ymis[i,k]] = tablamis[Ymis[i,k-1],Ymis[i,k]] +1
  }
}
tabla
tablamis

##################################################

##################################################
### HMM progression model
###
### Ordinal response 
### Continuos lanten process
### Non-decreasing process
###
### Run the  following lines of code
##################################################

##################################################

param <- c(
  "Beta" ,
  "Gama" ,
  "kappa" ,
  "Tau" ,
  "Sig2" ,
  "Sig" 
)

data <- list(
  Yobs = Ymis-1 ,
  X = X , 
  Z = Z ,
  N = N , 
  Xi = Xi ,
  K = K 
)

inits <- function(){	list(
  "Beta" = rnorm(2,0,0.01) , 
  "Gama" = rnorm(2,0,0.01) , 
  "Tau" = rep(1,2) , 
  "gam" = c(0,1,2) , 
  "Wobs" = Ymis-1.5  ,
  "Wtrue" = ifelse(Ymis>0,1,0)
)}
  
  
fit.ord <- jags.model("HMMprogressionOrdinal.bug", data, inits, n.chains=3)
  
update(fit.ord,2000)
sample.ord <- coda.samples(fit.ord, param, n.iter=2000, thin=2)
  
plot(sample.ord)
summary(sample.ord)
  
post.ord <- summary(sample.ord)
  
table1 <- cbind(post.ord$statistics[,1], post.ord$quantiles[,"50%"],
                post.ord$statistics[,2], post.ord$quantiles[,c("2.5%","97.5%")])
colnames(table1)[1:3] <- c("Mean", "Median", "SD") 
table1
#write.table(table1, file=paste0("table.txt"), quote=FALSE, sep=";", row.names=TRUE)

#typeof(sample.ord)
#save( sample.ord, file=paste0("sample.ord.RData"))
#load(file=paste0("sample.ord.RData"))
#summary(sample.ord)
  
  
##################################################

##################################################
##################################################
