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
### Aortic Aneurysm Progression Data
##################################################

##################################################
### R packages
### Load the R libraries
library(msm)
library(dplyr)
library(rjags)

##################################################
### ADDRESS
### Instructions: 
### Change the address where the data and codes are located. 
### setwd("HERE")
setwd("/Users/EstadisticaCiencias/Dropbox/HMM-Ordinal/HMMordinalCodes/")
getwd()

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
### Aortic Aneurysm Progression Data

data(aneur)
attach(aneur)
head(aneur)
help(aneur)

n_distinct(aneur$ptnum)
max(table(aneur$ptnum))
table(aneur$state)

N = 838   # subjects
K = 21   # times
J = 4   # categories

X = array(NA,dim=c(N,K,2)) ### (age,diam)
Y = array(NA,dim=c(N,K))

Ki = table(aneur$ptnum)
for(i in 1:N){
	aneur_i = aneur[aneur$ptnum==i,]
	for(k in 1:Ki[i]){
		Y[i,k] = aneur_i$state[k]
		X[i,k,1] = aneur_i$age[k]   
		X[i,k,2] = aneur_i$diam[k] 
	}  
} 
Xi = Y

##################################################
### Including all data

Y1 = Y 
X1 = X
N1 = N
Ki1 = Ki
Xi1 = Xi

##################################################
### Considering only data having more than one screen (state>1)

idx2 = c()
for(i in 1:N){
  if( sum(Y[i,1:Ki[i]])>Ki[i]){
    idx2 = c(idx2,i) 
  }
}

Y2 = Y[idx2,] 
X2 = X[idx2,,]
N2 = length(idx2)
Ki2 = Ki[idx2]
Xi2 = Xi[idx2,]

##################################################

##################################################
### Parameters, Initials and Data

param3.Mis <- c(
  "Beta" ,
  "kappa" ,
  "Tau" ,
  "Sig2" , 
  "Wtruenew" ,
  "Wobsnew" 
)

inits3.Mis <- function(){	list(
  "Beta" = rnorm(1,0,0.01) , 
  "Tau" = rep(1,J) , 
  "gam" = c(0,1,2) , 
  "Wobs" = Y2-1.5  ,
  "Wtrue" = ifelse(Y2>0,1,0)
)}

data3.Mis <- list(
  Yobs = Y2-1 ,
  X = X2[,,1] ,  
  N = N2 ,   
  Ki = Ki2 ,
  Xi = Xi2 ,
  Xnew = X2[c(67,80,119),,1] ,  
  Nnew = length(N2[c(67,80,119)]) ,   
  Kinew = Ki2[c(67,80,119)] ,
  Xinew = Xi2[c(67,80,119),] 
)

##################################################
### Fit the model

fit3.Mis <- jags.model("HMManeur.bug", data3.Mis, inits3.Mis, n.chains=1)

update(fit3.Mis,50000)
sample3.Mis <- coda.samples(fit3.Mis, param3.Mis, n.iter=50000, thin=5)


plot(sample3.Mis)
summary(sample3.Mis)
post3.Mis <- summary(sample3.Mis)

#typeof(sample3.Mis)
#save( sample3.Mis, file=paste0("Aneur.sample3.Mis.RData"))
#load(file=paste0("Aneur.sample3.Mis.RData"))
#summary(sample3.Mis)

##################################################
### Graphics

plot(0, xlim=c(60,79.5), ylim=c(10,63), type="n", main="Observed",xlab="Age at examination",ylab="Aortic diameter (mm)")

abline(h=c(30,45,55),lty=2,lwd=2)
text(61.5,25,"Stage 1")
text(61.5,37,"Stage 2")
text(61.5,50,"Stage 3")
text(61.5,59,"Stage 4")

points(X2[67,,1],X2[67,,2], pch=19,cex=0.8,col="red")
lines(X2[67,,1],X2[67,,2], lwd=1.5,col="red")
points(X2[80,,1],X2[80,,2], pch=17,cex=0.8,col="magenta")
lines(X2[80,,1],X2[80,,2], lwd=1.5,col="magenta")
points(X2[119,,1],X2[119,,2], pch=15,cex=0.8,col="blue")
lines(X2[119,,1],X2[119,,2], lwd=1.5,col="blue")

Ki2[c(67,80,119)]

legend("bottom",c("Subject 690","Subject 705","Subject 746"),pch=c(19,17,15),col=c("red","magenta","blue"),lty=1,ncol=3,cex=0.8)

##################################################
### Predictions

Ytrue67 = post3.Mis$statistics[c("Wtruenew[1,1]","Wtruenew[1,2]","Wtruenew[1,3]","Wtruenew[1,4]","Wtruenew[1,5]","Wtruenew[1,6]","Wtruenew[1,7]","Wtruenew[1,8]"), 1]
Ytrue80 = post3.Mis$statistics[c("Wtruenew[2,1]","Wtruenew[2,2]","Wtruenew[2,3]","Wtruenew[2,4]","Wtruenew[2,5]","Wtruenew[2,6]","Wtruenew[2,7]","Wtruenew[2,8]","Wtruenew[2,9]","Wtruenew[2,10]","Wtruenew[2,11]","Wtruenew[2,12]","Wtruenew[2,13]","Wtruenew[2,14]","Wtruenew[2,15]"), 1]
Ytrue119 = post3.Mis$statistics[c("Wtruenew[3,1]","Wtruenew[3,2]","Wtruenew[3,3]","Wtruenew[3,4]","Wtruenew[3,5]","Wtruenew[3,6]","Wtruenew[3,7]","Wtruenew[3,8]","Wtruenew[3,9]","Wtruenew[3,10]","Wtruenew[3,11]","Wtruenew[3,12]","Wtruenew[3,13]","Wtruenew[3,14]","Wtruenew[3,15]"), 1]

kappa = post3.Mis$statistics[c("kappa[1]","kappa[2]","kappa[3]"), 1]


plot(0, xlim=c(60,79.5), ylim=c(-1.2,2.7), type="n", main="Estimated", xlab="Age at examination",ylab="Stages", labels=FALSE,tick=FALSE)
axis(1,at=c(60,65,70,75,80),lab=c(60,65,70,75,80))

abline(h=kappa,lty=2,lwd=2)
text(61.5,-0.4,"Stage 1")
text(61.5,0.6,"Stage 2")
text(61.5,1.7,"Stage 3")
text(61.5,2.4,"Stage 4")

points(X2[67,1:8,1],Ytrue67, pch=19,cex=0.8,col="red")
lines(X2[67,1:8,1],Ytrue67, lwd=1.5,col="red")
points(X2[80,1:15,1],Ytrue80, pch=17,cex=0.8,col="magenta")
lines(X2[80,1:15,1],Ytrue80, lwd=1.5,col="magenta")
points(X2[119,1:15,1],Ytrue80, pch=15,cex=0.8,col="blue")
lines(X2[119,1:15,1],Ytrue80, lwd=1.5,col="blue")

legend("bottom",c("Subject 690","Subject 705","Subject 746"),pch=c(19,17,15),col=c("red","magenta","blue"),lty=1,ncol=3,cex=0.8)

##################################################

##################################################
##################################################

