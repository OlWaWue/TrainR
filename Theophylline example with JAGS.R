

library(R2jags)
library(ggplot2)
data("Theoph")

Theoph$Dose <- Theoph$Dose*Theoph$Wt

colnames(Theoph) <- c("ID", "WT", "AMT", "TIME", "DV")

ggplot(data=Theoph) + geom_point(aes(x=TIME, y=DV, colour=ID ))

nobs <- nrow(Theoph)
start <- (1:nobs)[!duplicated(Theoph$ID)]
end = c(start[-1] - 1, nobs)
nsub = length(unique(Theoph$ID))


bugsdata <- list(
  nobs = nobs,
  nsub = nsub,
  start = start,
  end = end,
  weight = Theoph$WT,
  time = Theoph$TIME,
  amt = Theoph$AMT,
  logCobs = ifelse(Theoph$DV <= 0, NA, log(Theoph$DV)),
  omegaInvPrior = 3 * diag(rep(0.05, 3))
)

model <- function() {
    
    for(i in 1:nsub){
      
      ## Inter-individual variability
      
      ## Individual parameters, i.e., parameters conditioned on data for the ith individual
      logtheta[i, 1:3] ~ dmnorm(logthetaMean[i, 1:3], omegaInv[1:3, 1:3])
      logthetaMean[i, 1] <- logkeHat# + 0.75*log(weight[start[i]]/70) # ke
      logthetaMean[i, 2] <- logV1Hat# + log(weight[start[i]]/70) # V1
      logthetaMean[i, 3] <- log(kaHat)     # ka
      
      theta[i,4] <- 1                     # F1

      
      ## Simulated parameters for hypothetical new individuals
      logthetaPred[i, 1:3] ~ dmnorm(logthetaMean[i, 1:3], omegaInv[1:3, 1:3])
      thetaPred[i,4] <- 1                 # F1

      
      for(j in 1:3){
        log(theta[i,j]) <- logtheta[i,j]
        log(thetaPred[i,j]) <- logthetaPred[i,j]
      }
      
      ## Calculate amounts in each compartment for all event times for the ith individual
      
      
      
      for (k in start[i]:end[i]) {
        
      ## Individual predictions
      xhat[k] <- theta[i,4]*amt[k]/theta[i,2] *(theta[i,3]/(theta[i,3]- theta[i,1]))*( (exp(-theta[i,1]*(time[k]) )) - (exp(-theta[i,3]*(time[k])))  )

      ## Population predictions
      xhatPred[k] <- thetaPred[i,4]*amt[k]/thetaPred[i,2] *(thetaPred[i,3]/(thetaPred[i,3]-thetaPred[i,1]))*( (exp(-thetaPred[i,1]*(time[k]) )) - (exp(-thetaPred[i,3]*(time[k])))  )
      }
    }
    
    for(i in 1:nobs){
      
      logCobs[i] ~ dnorm(logCHat[i], tau)# Likelihood
      logCobsCond[i] ~ dnorm(logCHat[i], tau) # Individual prediction
      CHat[i] <- xhat[i] 
      logCHat[i] <- log(max(CHat[i], eps))
      
      logCobsPred[i] ~ dnorm(logCHatPred[i], tau) # Population prediction
      CHatPred[i] <- xhatPred[i] 
      logCHatPred[i] <- log(max(CHatPred[i], eps))
      
    }
    
    ## Prior distributions
    
    logkeHat ~ dnorm(0, 1.0E-6)
    logkaHat ~ dnorm(0, 1.0E-6)
    logV1Hat ~ dnorm(0, 1.0E-6)
    log(keHat) <- logkeHat
    log(kaHat) <- logkaHat
    log(V1Hat) <- logV1Hat
    
    ##  tau <- 1 / (sigma * sigma)
    ##  sigma ~ dunif(0, 1000)
    tau ~ dgamma(0.001, 0.001)
    sigma <- sqrt(1/tau)
    omegaInv[1:3, 1:3] ~ dwish(omegaInvPrior[1:3, 1:3], 3)
    omega[1:3, 1:3] <- inverse(omegaInv[1:3, 1:3])
    eps <- 1.0E-6
    
}

## Specify the variables for which you want history and density plots
parametersToPlot <- c("keHat", "kaHat", "V1Hat",
                      "sigma", "omega", "CHat")

## Additional variables to monitor
otherRVs <- c("logCobsCond", "logCobsPred", "theta")

parameters <- c(parametersToPlot,otherRVs)

## create initial estimates
bugsinit <- function() {
  list(logkeHat = rnorm(1, log(2), 0.2),
       logkaHat = rnorm(1, log(2), 0.2),
       logV1Hat = rnorm(1, log(70), 0.2),
       omegaInv = solve(diag(exp(2 * rnorm(3, log(0.25), 0.5)))),
       tau = 1/(runif(1, 0.1, 2)^2)
       ##         sigma = runif(1, 0.1, 2)
  )
}

jagsfit <- jags(bugsdata, model=model, inits = bugsinit,
                 parameters=parameters,progress.bar="gui",
                n.chains=4,n.iter=1400,n.burnin=500,n.thin=2)



pk_mod <- function(time=1:24, amt=100, par=c(ka=2,ke=0.1,V=2)){
  conc <- amt/par[["V"]]*(par[["ka"]]/(par[["ka"]]-par[["ke"]]))*(exp(-par[["ke"]]*time)-exp(-par[["ka"]]*time))
  
  return(conc)
}

pop_V <- jagsfit$BUGSoutput$mean$V1Hat
pop_ke <- jagsfit$BUGSoutput$mean$keHat
pop_ka <- jagsfit$BUGSoutput$mean$kaHat



Theoph$POP_PRED <- pk_mod(Theoph$TIME, mean(Theoph$AMT), par=c(ka=ka,
                                                               ke=ke,
                                                               V=V))


ggplot(data=Theoph) + geom_point(aes(x=TIME, y=DV)) + 
  geom_line(aes(x=TIME, y=POP_PRED), colour="red")  + facet_wrap(.~ID)


