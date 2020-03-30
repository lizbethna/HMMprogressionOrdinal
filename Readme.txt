##################################################

Paper: 
A hidden Markov model to address measurement errors in ordinal response scale and non-decreasing process


Authors:
Lizbeth Naranjo (1), Luz Judith R. Esparza (2), Carlos J. Perez (3).

(1) Departamento de Matemáticas, Facultad de Ciencias, Universidad Nacional Autonoma de Mexico (UNAM), Mexico

(2) Catedra CONACyT - Universidad Autonoma de Aguascalientes, Mexico

(3) Departamento de Matematicas, Facultad de Veterinaria, Universidad de Extremadura, Spain


Journal: Mathematics
Submitted. Under Revision. 


##################################################

Instructions to run the codes in R and JAGS are provided. 
The codes are applied to obtain a similar analysis as in Section 4 ‘Simulation example’, but without cross-validation, and Section 5 ‘Aortic aneurysm progression’. 

##################################################

##################################################
FILES 

For Section 4 ‘Simulation example’
The file ‘HMMprogressionOrdinal.R’ contains the R code. The JAGS code is run from this R file.

The file ‘HMMprogressionOrdinal.bug' contains the JAGS model. 

For Section 5 ‘Aortic aneurysm progression’
The file ‘HMManeur.R’ contains the R code. The JAGS code is run from this R file.. 

The file ‘HMManeur.bug' contains the JAGS model. 


##################################################

To run the files, do the following.
 
1.- Download JAGS from www.mcmc-jags.sourceforge.net/

2.- Install the packages necessary to run the R file. 
These are indicated in the R file. 

3.- Change the address indicated in ‘setwd()’. 
setwd("HERE"). 
This is the address where the file ‘HMMprogressionOrdinal.bug’ and ‘HMManeur.bug' are in.

##################################################


