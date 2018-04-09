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

###################################################################################

 model{ 

# parameter model- uninformed priors

  k~dbeta(1,1)
  sigmap~dgamma(0.001,0.001)
      
# likelihood 

  for( i in 1:n) {
     
      Mt[i] ~ dgamma(mu[i]^2/sigmap^2,mu[i]/sigmap^2) # moment matching 
      mu[i]<- Mo[i]*exp(-k*t[i])  # process model
      Mt.new[i] ~ dgamma(mu[i]^2/sigmap^2,mu[i]/sigmap^2)   
        
} # end of loop



  # Bayesian p values- pvalue.cv is very large (close to 1) or very small (close to 0) model fails to
   #adequately rep distribution, extreme would be <0.10 or >0.90, good representation bayesian p value #near 0.50

cv.y <- sd(Mt[ ])/mean(Mt[ ])
cv.y.new <- sd(Mt.new[])/mean(Mt.new[ ])
pvalue.cv <- step(cv.y.new-cv.y)

mean.y <-mean(Mt[])
mean.y.new <-mean(Mt.new[])
pvalue.mean <-step(mean.y.new - mean.y)


for (j in 1:n){
  sq[j]<-(Mt[j]-mu[j])^2
  sq.new[j]<-(Mt.new[j]-mu[j])^2
  
}

fit <- sum(sq[])
fit.new <- sum(sq.new[])
pvalue <- step(fit-fit.new)     # bayesan p value


}

 
 