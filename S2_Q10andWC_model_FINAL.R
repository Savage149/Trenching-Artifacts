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


###############################################################################################
## ALTERING CODE FOR DIFFERENT TEMPERATURE AND MOISTURE MODELS
## LINE 29 CONTAINS THE MODEL CODE AND WOULD NEED TO BE ADJUSTED TO THE NEW MODEL
## UNDER PARAMETER MODEL; THIS IS WHERE THE USER MUST DEFINE NEW PROBABILITY DISTRIBUTIONS 
## FOR THIER SPECIFIC MODEL (LINES 16-18)

## ONCE THAT IS CHANGED THE USER MUST CHANGE THE PARAMETERS IN S2_Q10andWC_jags_FINAL.R
## IN S2_Q10andWC_jag_FINAL
## LINE 20-23 INITS MUST BE CHANGED TO NEW PARAMETERS AND USER DETERMINED START POINTS
## LINE 27- MAKE SURE REFERENCES TO CORRECT RS, SOILT AND SOILMOISTURE FROM INPUT FILE MYDATA
## LINE 47-TO END- MAKE SURE TO UPDATE MODEL VARIABLE NAMES- FOLLOW EXAMPLE FROM THIS CODE
## AND CHANGE PARAMETERS TO REFLECT USERS NEW MODEL

 model{ 

# parameter model
  Rref~dgamma(0.05,0.001)  # this is a mean Rref of 50 mg C m-2 hr-1,gamma distribution used becasuse Rref must be positive 
  Q~dgamma(0.002,0.001)   #this is mean Q prior of 2, gamma distribution used becasuse Q must be positive 
  B~dbeta(5,1)            #  this is a B mean of 0.83 , beta distribution used becasue B ranges between 0-1 
  tau~dgamma(0.001,0.001)
  sigma<-1/sqrt(tau)
    
  # know priors for Rref, Q and B based on work from Savage et al. 2013 and Savage et al. 2009 

  
  for( i in 1:n) {
    
      RS[i] ~ dnorm(log(mu[i]),tau) # log transform RS data prior to use normal distribution     
      mu[i]<- Rref*Q^((soilT[i]-10)/10)*(B^VSM[i]) # process model
      RS.new[i] ~ dnorm(log(mu[i]),tau)  
      sq[i]<-(RS[i]-log(mu[i]))^2
      sq.new[i]<-(RS.new[i]-log(mu[i]))^2
      
        
       } # end of loop

   # predictive check
  # Bayesian p values- pvalue.cv is very large (close to 1) or very small (close to 0) model fails to 
  #adequately rep distribution, extreme would be <0.10 or >0.90, good representation bayesian p value #near 0.50
 
cv.y <- sd(RS[ ])/mean(RS[ ])
cv.y.new <- sd(RS.new[])/mean(RS.new[ ])
pvalue.cv <- step(cv.y.new-cv.y)
  
mean.y <-mean(RS[])
mean.y.new <-mean(RS.new[])
pvalue.mean <-step(mean.y.new - mean.y)
  

fit <- sum(sq[])
fit.new <- sum(sq.new[])
pvalue <- step(fit-fit.new)     # Bayesan p value

 }
 
 # end of model
