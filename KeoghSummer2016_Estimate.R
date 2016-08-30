############################
# KeoghSummer2016_Estimate.R
# Keogh2016.Rproj

# Calculates a 2 or 3-pass depletion estimate for 2016 Keogh data

# Files sourced:
# Files produced:

# Created August 26, 2016
# A Putt
#############################

# Methods are from Ogle FSA Package
library(FSA)

# Source the data
source("KeoghSummer2016_DataUpload.R")

#########################################
# 2-Pass SeberLeCren ####################
#########################################

# The Method Matt referenced in Seber and Cren
# 2-pass removal assuming a large proportion of the population is removed
# Modified Leslie method, modified by Seber and le Cren 1967 

SeberLeCren <- function(countvec) {
  C1 <- countvec[1]
  C2 <- countvec[2]
  q <- C2/C1 
  p <- 1-q # Probability any one fish is caught
  N <- (C1*C1)/(C1-C2) # Population estimate
  StErrN <- sqrt( ( C1^2*C2^2*(C1+C2) ) / ( (C1-C2)^4 ) )
  upperCI <- N + StErrN*1.96
  lowerCI <- N - StErrN*1.96
  result <- data.frame(q=q,p=p,N=N,StErrN=StErrN,upperCI=upperCI,lowerCI=lowerCI)
  return(result)
}

########################################
# 3-Pass K-Pass Carle Strub Model ######
########################################
# 3-pass removal with costant effort per pass

# 1. Zippin Method
# This is the Zippin method, which doesn't consider Bayesian priors
# Equation 3 in Carle Strub

KPass <- function(C1,C2,C3) {
  ct <- c(C1,C2,C3)
  k <- length(ct) # Number of removals (in this case it will be three as we are only looking at 3-pass)
  T <- sum(ct) # Total catch
  i <- seq(1,k) # Vector of pass numbers
  X <- sum((k-i)*ct) # I can't figure out what this actually is
  mle <- function(N) { (N+0.5)*((k*N-X-T)^k) - (N-T+0.5)*((k*N-X)^k) } # The first integer that results in a negative function result is the population size
  NVector <- seq(T,500) # I'm not sure how to solve the above rule iteratively so I'll do it the long way
  mleResult <- data.frame(N=NVector,mle=mle(NVector)) # The first negative value is the population estimate
  FirstNegative <- min(which(mleResult$mle <= 0)) # Find the negative rows then take the minimum
  N <- NVector[FirstNegative]
}

# 2. Carle Strub Method
# The Method Matt referenced in Carle and Strub 1978 (Equation 7)
# This can also be done using the FSA package by Ogle
# The FSA package is faster than my above function and also easily computes error and 95% CI
# requires a vector of catches as the first argument
# If we set alpha and beta to 1 we just get an uninformative prior
# Could also be used for 2-pass

# removal(catch vector,method="CarleStrub",alpha=1,beta=1)
# removal(catch vector,method="Zippin")

###################################
# Calculate Estimates #############
###################################

head(counts)

# Because of the table format and short length of the table it is easiest to run a loop

# Create a vector of unique sites
site <- levels(counts$uniquesite)

# Create an empty table for the loop to populate
estimateList <- list()

for (i in 1:length(site)) {
  OneSite <- subset(counts,uniquesite==site[i]) # Pull one site
  if (nrow(OneSite)==2) { # Run the Seber Le Cren Method for 2-pass
    cof <- SeberLeCren(OneSite$cof.count)
    cop <- SeberLeCren(OneSite$cop.count)
    shf <- SeberLeCren(OneSite$shf.count)
    shp <- SeberLeCren(OneSite$shp.count)
    Estimate <- c(cof$N,cop$N,shf$N,shp$N)
    StErr    <- c(cof$StErrN,cop$StErrN,shf$StErrN,shp$StErrN)
    upperCI  <- c(cof$upperCI,cop$upperCI,shf$upperCI,shp$upperCI)
    lowerCI  <- c(cof$lowerCI,cop$lowerCI,shf$lowerCI,shp$lowerCI)
    Method   <- rep("2-Pass Seber Le Cren",4)
    
  } else if (nrow(OneSite)==3) { # Run the Carle Strub method for 3-pass
    cof <- removal(OneSite$cof.count,method="CarleStrub",alpha=1,beta=1)$est
    cop <- removal(OneSite$cop.count,method="CarleStrub",alpha=1,beta=1)$est
    shf <- removal(OneSite$shf.count,method="CarleStrub",alpha=1,beta=1)$est
    shp <- removal(OneSite$shp.count,method="CarleStrub",alpha=1,beta=1)$est
    Estimate <- as.vector(c(cof[1],cop[1],shf[1],shp[1]))
    StErr    <- as.vector(c(cof[2],cop[2],shf[2],shp[2]))
    lowerCI  <- as.vector(c(cof[3],cop[3],shf[3],shp[3]))
    upperCI  <- as.vector(c(cof[4],cop[4],shf[4],shp[4]))
    Method   <- rep("3-Pass Carle Strub",4)
  }
  
  result <- data.frame(Site=rep(site[i],4),Species=c("cof","cop","shf","shp"),Estimate=Estimate,
                       StErr=StErr,lowerCI=lowerCI,upperCI=upperCI,Method=Method)
  estimateList[[i]] <- result
}
  

# Turn the list into a data.frame
EstimatesDF <- do.call("rbind",estimateList)

# Pair Estimates with site areas and calculate densities
areas <- data.frame(Site=locations$uniquesite,Area=locations$Area_m2) # Create data frame to merge
Estimates <- merge(EstimatesDF,areas,by="Site") # Merge the two frames by unique site
Estimates$density_per_m2 <- Estimates$Estimate/Estimates$Area # Caluculate density/m2
  
write.csv(Estimates,"R_ResultOutput/Estimates.csv",row.names=FALSE)  
  