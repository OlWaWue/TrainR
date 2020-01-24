library("R2jags")

a = 10
b = 20
slope = a/b

numberOfRealdata = 100

xsamples <- seq(length=numberOfRealdata, from=1, to=numberOfRealdata)

theoretical_hockeyStick <- ifelse(xsamples < a, xsamples * slope, b)

noisy_hockeyStickSamples <- theoretical_hockeyStick+runif(numberOfRealdata, -1, +1)

SD.slope = 0.01
SD.a = 0.01
SD.b = 0.01
N= numberOfRealdata

jags.data <- list("N", "a", "b", "slope", "SD.a", "SD.b", "SD.slope",
                  "xsamples", "noisy_hockeyStickSamples")

jags.parms <- c("random_a", "random_b", "random_slope")

Model = "
model{

random_slopetau <- pow(SD.slope, -2)
random_a_tau <- pow(SD.a, -2)
random_b_tau <- pow(SD.b, -2)

random_slope ~ dnorm(slope, random_slopetau)
random_a ~ dnorm(a, random_a_tau)
random_b ~ dnorm(b, random_b_tau)

  for (j in 1:N){
    y[j] <- ifelse(xsamples[j] < random_a, xsamples[j] * random_slope, random_b)
    noisy_hockeyStickSamples[j] ~ dnorm(y[j], random_b_tau)
  }

}"

library("coda")

JAGSFILE="r2ssb.bug"

cat(Model, file = JAGSFILE)

Nchains = 2
Nburnin = 100
Niter = 1000
Nthin = 10

jagsfit <- jags(data=jags.data,
                working.directory = NULL,
                inits = NULL,
                jags.parms,
                model.file = JAGSFILE,
                n.chains=Nchains,
                n.thin = Nthin,
                n.iter = Niter,
                n.burnin = Nburnin)

random_a_samples <- jagsfit$BUGSoutput$sims.list$random_a
random_b_samples <- jagsfit$BUGSoutput$sims.list$random_b
random_slope <- jagsfit$BUGSoutput$sims.list$random_slope

a_best <- mean(random_a_samples)
SD.a_best <- apply(as.matrix(a_best), 2, sd)
b_best <- mean(random_b_samples)
SD.b_best <- apply(as.matrix(b_best), 2, sd)

cat("Best estimate for a:", a_best, "\n")
cat("Best estimate for b:", b_best, "\n")

plot(xsamples, noisy_hockeyStickSamples, col="blue")
abline(v=a, lty=3, col="red", lwd=1.5)
text(x=a+10, y=9, "True value of a")
lines(x=c(0, a_best, numberOfRealdata), y=c(0, b_best, b_best), col="green")
text(x=a+40, y=17, "Estimated Hockey-Stick function")
