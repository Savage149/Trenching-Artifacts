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

#######################################################################################
# read in Fluxdata_example.csv
# this data is a subset of the data used for the manuscript and is intended as an example
# makes sure Q10andWC_model_FINAL.R is in the same directory as this R code
##IMPORTANT- LATEST VERSION OF RSTUDIO Version 1.0.153 DO NOT RUN RJAGS CORRECTLY, YOU WILL
### GET A CODA ERROR. YOU CAN RUN RJAGS WITH OLDER VERSION OR RUN FROM R COMMAND LINE ONLY.  
### HOPEFULLY THIS WILL GET FIXED SOON, CONTINUE TO CHECK UPDATES.

rm(list=ls())
mydata=read.table(file.choose(),sep=",", header=T)

install.packages('rjags')
library(rjags)
library(lattice)

##########################################################################################


log_Flux<-log(mydata$Flux)

# initial parameters- note they must always be a list of lists
# this is a list with 3 chains, given differing starting points

inits= list(
  list(Rref=50, Q=2.1,B=0.8,tau=0.01),
  list(Rref=150, Q=1.9, B=0.5, tau=0.001),
  list(Rref=100, Q=2.0, B=0.9, tau=0.001))
   

n.xy = nrow(mydata) 
data=list(n=n.xy, RS=log_Flux,soilT=mydata$soilT, VSM=mydata$VSM_optT)

# number of iterations

n.adapt=5000
n.update = 10000
n.iter = 10000

##call to JAGS- this sets up MCMC sampling

jm=jags.model("Q10andWC_model_FINAL.R",data=data, inits,n.chains=length(inits), n.adapt=n.adapt)

#burnin chain- these data are not used for posterior estimates
update(jm, n.iter=n.update)

# this module needs to be loaded before coda and jags statement
load.module("dic")  

#generate coda object- the variables we want postierior distribution on

zm=coda.samples(jm,variable.names=c("Rref","Q","B", "sigma"), n.iter=n.iter, n.thin=1)
zj=jags.samples(jm,variable.names=c("Rref","Q","B","sigma","mu", "pvalue", "RS.new", "fit", 'fit.new'), n.iter=n.iter, n.thin=1)


# combine output into a dataframe and write MCMC output for further analysis- change FILENAME
df=as.data.frame(rbind(zm[[1]],zm[[2]],zm[[3]]))
write.table(df,"FILENAME",sep="\t",col.names=NA)

# SUMMARY PLOTS
acfplot(zm)
summary(zm)
plot(zm)

#######################################################################################################
# checking for convergence

q.coda = as.mcmc.list(zj$Q)
b.coda = as.mcmc.list(zj$B)
R.coda = as.mcmc.list(zj$Rref)
plot(q.coda)
xyplot(q.coda)
densityplot(q.coda)

# Gelman and Rubin diagnostics  - determines when the chain has forgotten thier inital values

gelman.diag(zm)     # convergence of point.est approaches 1- more iterations if 95% CI >1.05


#####POSTERIOR PREDICTIVE CHECK ####################################
### values range between 0-1, we want mid range values, not extremes <0.1 or >0.9, this indicates bad model choice

zj$pvalue  # Bayesian p value



