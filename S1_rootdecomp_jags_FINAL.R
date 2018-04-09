### Biogeochemistry
### Partitioning Soil Respiration: Quantifying the Artifacts of the Trenching Method
###Authors: K E Savage1, E A Davidson2, R Z Abramoff3,4, A C Finzi3, M-A Giasson3
### Affiliations
### 1 Woods Hole Research Center. 149 Woods Hole Rd., Falmouth, MA, 02540
### 2 University of Maryland Center for Environmental Science. 301 Braddock Road, Frostburg, MD 21532
### 3Boston University. 5 Cummington Mall, Boston, MA 02215
### 4Lawrence Berkley National Laboratory. One Cyclotron Road, MS74R316C, Berkeley, CA 94720

###Corresponding author
###Kathleen Savage
###149 Woods Hole Rd
###Falmouth, MA
###02540
###508-444-1542
###savage@whrc.org

####################################################################################


# read in Rootdecomp.csv
# this data is a subset of the data used for the manuscript and is intended as an example
# makes sure decomp_model_FINAL.R is in the same directory as this R code
##IMPORTANT- LATEST VERSION OF RSTUDIO Version 1.0.153 DO NOT RUN RJAGS CORRECTLY, YOU WILL
### GET A CODA ERROR. YOU CAN RUN RJAGS WITH OLDER VERSION OR RUN FROM R COMMAND LINE ONLY.  
### HOPEFULLY THIS WILL GET FIXED SOON, CONTINUE TO CHECK UPDATES.

rm(list=ls())
mydata=read.table(file.choose(),sep=",", header=T)

###########################################################################################

# initial parameters- note they must always be a list of lists
install.packages('rjags')
library(rjags)
library(lattice)

# initial parameters- note they must always be a list of lists
# this is a list with 3 chains, given differing starting points

inits= list(
  list(k=0.0007,sigmap=0.01),
  list(k=0.002, sigmap=0.001),
  list(k=0.1,sigmap=0.0001)) 


n.xy = nrow(mydata) 
data=list(n=n.xy, Mt=mydata$Mt,Mo=mydata$Mo, t=mydata$t)

# number of iterations

n.adapt=5000
n.update = 10000
n.iter = 10000

##call to JAGS- this sets up MCMC sampling

jm=jags.model("decomp_model_FINAL.R",data=data, inits,n.chains=length(inits), n.adapt=n.adapt)

#burnin chain

update(jm, n.iter=n.update)

# this module needs to be loaded before coda and jags statement

load.module("dic")  

#generate coda object- the variables we want postierior distribution on
# deviance- sum of the deviance of all observed stochastic nodes
# coda- these are the variables we want postieror distribution on 

zm=coda.samples(jm,variable.names=c("k","sigmap"), n.iter=n.iter, n.thin=1)
zj=jags.samples(jm,variable.names=c("k","sigmap","mu","pvalue"), n.iter=n.iter, n.thin=1)

# combine output into a dataframe and write MCMC output for further analysis- change FILENAME
df=as.data.frame(rbind(zm[[1]],zm[[2]]))
write.table(df,"FILENAME",sep="\t",col.names=NA)

# SUMMARY PLOTS
acfplot(zm)
summary(zm)
plot(zm)

#######################################################################################################
# checking for convergence- convert jags to coda object for this

k.coda = as.mcmc.list(zj$k)
plot(k.coda)
xyplot(k.coda)
densityplot(k.coda)

# Gelman and Rubin diagnostics  - determines when the chain has forgotten thier inital values

gelman.diag(zm)       # convergence of point.est approaches 1- more iterations if 95% CI >1.05
                          # note gelman.diag(zm) gives all monitored parameters  




