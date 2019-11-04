IV <- data.frame(
  time=c(0.33,0.5,0.67,1.5,2,4,6,10,16,24,32,48),
  conc=c(14.7,12.6,11,9,8.2,7.9,6.6,6.2,4.6,3.2,2.3,1.2))
Oral <- data.frame(
  time=c(0.5,1,1.5,2,4,6,10,16,24,32,48),
  conc=c(2.4,3.8,4.2,4.6,8.1,5.8,5.1,4.1,3,2.3,1.3))

plot(conc~time,data=IV[1:12,],log='y',col='blue',type='l')
with(Oral,lines(y=conc,x=time,col='red'))

library(R2jags)

datain <- list(
  timeIV=IV$time,
  concIV=IV$conc,
  nIV=nrow(IV),
  doseIV=500,
  timeO=Oral$time,
  concO=Oral$conc,
  nO=nrow(Oral),
  doseO=500)

model1 <- function() {
  tau <- 1/pow(sigma,2)
  sigma ~ dunif(0,100)
  # IV part
  kIV ~dnorm(.4,1)
  cIV ~ dlnorm(1,.01)
  for (i in 1:nIV) {
    predIV[i] <- c0*exp(-k*timeIV[i]) +cIV*exp(-kIV*timeIV[i])
    concIV[i] ~ dnorm(predIV[i],tau)
  }
  c0 <- doseIV/V
  V ~ dlnorm(2,.01)
  k <- CL/V
  CL ~ dlnorm(1,.01)
  AUCIV <- doseIV/CL+cIV/kIV
  # oral part
  for (i in 1:nO) {
    predO[i] <- c0star*(exp(-k*timeO[i])-exp(-ka*timeO[i]))
    concO[i] ~ dnorm(predO[i],tau)
  }
  c0star <- doseO*(ka/(ka-k))*F/V
  AUCO <- c0star/k
  F ~ dunif(0,1)
  ka ~dnorm(.4,1)
  ta0_5 <- log(2) /ka
  t0_5 <- log(2)/k
}

parameters <- c('k','AUCIV','CL','V','t0_5','c0',
                'AUCO','F','ta0_5','ka','c0star',
                'kIV','cIV','predIV','predO')
inits <- function() 
  list(
    sigma=rnorm(1,1,.1),
    V=rnorm(1,25,1),
    kIV=rnorm(1,1,.1),
    cIV = rnorm(1,10,1),
    CL = rnorm(1,5,.5),
    F = runif(1,0.8,1),
    ka = rnorm(1,1,.1)
  )
jagsfit <- jags(datain, model=model1, 
                inits=inits, parameters=parameters,progress.bar="gui",
                n.chains=4,n.iter=14000,n.burnin=5000,n.thin=2)


jagsfit$BUGSoutput

cin <- c(IV$conc,Oral$conc)
mIV <- sapply(1:ncol(jagsfit$BUGSoutput$sims.list$predIV),
              function(x) {
                mean(jagsfit$BUGSoutput$sims.list$predIV[,x])
              }
)
mO  <- sapply(1:ncol(jagsfit$BUGSoutput$sims.list$predO),
              function(x) {
                mean(jagsfit$BUGSoutput$sims.list$predO[,x])
              }
)      
plot(x=c(IV$conc,Oral$conc),y=c(IV$conc,Oral$conc)-c(mIV,mO))


tpred <- 0:48
simO <- sapply(1:jagsfit$BUGSoutput$n.sims,function(x) {
  jagsfit$BUGSoutput$sims.list$c0star[x]*
    (exp(-jagsfit$BUGSoutput$sims.list$k[x]*tpred)
     -exp(-jagsfit$BUGSoutput$sims.list$ka[x]*tpred))
})
predO <- data.frame(mean=rowMeans(simO),sd=apply(simO,1,sd))

simIV <- sapply(1:jagsfit$BUGSoutput$n.sims,function(x) {
  jagsfit$BUGSoutput$sims.list$c0[x]*
    (exp(-jagsfit$BUGSoutput$sims.list$k[x]*tpred)
     +exp(-jagsfit$BUGSoutput$sims.list$kIV[x]*tpred))
})
predIV <- data.frame(mean=rowMeans(simIV),sd=apply(simIV,1,sd))

plot(conc~time,data=IV[1:12,],log='y',col='blue',type='p',ylab='Concentration')
with(Oral,points(y=conc,x=time,col='red'))
lines(y=predO$mean,x=tpred,col='red',lty=1)
lines(y=predO$mean+1.96*predO$sd,x=tpred,col='red',lty=2)
lines(y=predO$mean-1.96*predO$sd,x=tpred,col='red',lty=2)
lines(y=predIV$mean,x=tpred,col='blue',lty=1)
lines(y=predIV$mean+1.96*predIV$sd,x=tpred,col='blue',lty=2)
lines(y=predIV$mean-1.96*predIV$sd,x=tpred,col='blue',lty=2)



###
# With STAN
###


library('rstan');

source('data_stan.R');

fit <- stan('iv_oral.stan', 
            data=c("nIV","nOral","doseIV","doseOral", 
                   "timeIV","concIV","timeOral","concOral"),
            chains=4, warmup=5000, iter=14000)


