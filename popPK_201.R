library(mrgsolve)
library(tidyr)
library(dplyr)
library(ggplot2)

mod <- mread("Examplotinib")

set.seed(16111983)

patient_data <- data.frame(ID=1:20, SEX=rbinom(20, 1, 0.5), WT=rnorm(20, 75, 15), CLCR=rnorm(20,115, 35))



res <- mod %>% ev(ID=1:20, amt=1000, addl=0, ii=24) %>% idata_set(patient_data) %>% carry_out(amt,evid, addl, ii) %>% mrgsim(end=24)

plot(res)

dat <- as.data.frame(res)

ggplot(dat) + geom_point(aes(x=time, y=DV, colour=as.factor(SEX) )) + facet_wrap(dat$ID)


sample_dat <- filter(dat, dat$time %in% c(0,0.25, 0.5, 0.75, 1, 1.25, 2, 2.5, 4, 6, 8, 12, 24))



ggplot(sample_dat) + geom_point(aes(x=time, y=DV)) + facet_wrap(sample_dat$ID)

library(writexl) 

write_xlsx(sample_dat, path="./pat_dat.xlsx")

library(readxl)

data <- read_excel("./pat_dat.xlsx")


data$cmt <- 2

data$mdv <- 0

data$cmt <- ifelse(data$evid == 1, 1, 2)

data$mdv <- ifelse(data$evid == 1, 1, 0)

data$LNDV <- log(data$DV)



data <- data[-which(data$time==0 & data$evid==0),]

write.csv(data, "pat_dat.csv", row.names = F)

data <- read.csv("pat_dat.csv")

d <- nmDataConvert(data);





my2CptModel <- function() {
  ini({
    tka <- log(0.22)
    tcl <- log(3.395)
    tv2  <- log(35)
    tv3  <- log(450)
    tq   <- log(5)
    tlagtime <- log(0.5)
    eta.ka ~ 1
    eta.cl ~ 1
    eta.v2 ~ 1
    eta.v3 ~ 1
    eta.q ~ 1
    eta.tlag ~1
    add.err <- 0.075
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v2 <- exp(tv2 + eta.v2)
    v3 <- exp(tv3 + eta.v3)
    q  <- exp(tq + eta.q)
    tlag <- exp(tlagtime+eta.tlag)
    d/dt(depot) = -ka * depot
    d/dt(center) = ka * depot - cl / v2 * center + q/v3 * periph - q/v2 * center
    d/dt(periph) = q/v2 * center - q/v3 * periph
    cp = center / v2
    cp ~ add(add.err)
  })
}


library(nlmixr)
library(ggplot2)



ggplot(d) + geom_point(aes(x=TIME, y=DV)) + facet_wrap(d$ID)

ggplot(d) + geom_point(aes(x=TIME, y=LNDV)) + facet_wrap(d$ID)

my2CptFit <- nlmixr(my2CptModel, d, est="saem", control=saemControl(nBurn = 200, nEm=300))

plot(my2CptFit)

ggplot(my2CptFit) + geom_point(aes(x=TIME, y=DV, colour=ID)) + 
  geom_line(aes(x=TIME, y=IPRED, colour=ID)) + facet_wrap(my2CptFit$ID)

xpdb <- xpose.nlmixr::xpose_data_nlmixr(my2CptFit)

library(xpose.nlmixr)
library(xpose)

dv_vs_ipred(xpdb)

xpose::eta_distrib(xpdb)



library(nlmixr)
library(ggplot2)
library(xpose.nlmixr)
#-----------------------------------------------------------
#Exploratory plots

df = Theoph

ggplot(df, aes(x=Time, y=conc, color=Subject)) +
  geom_point() +
  geom_line()

ggplot(df, aes(x=Time, y=conc)) +
  geom_point() +
  stat_smooth() +
  scale_y_log10()

#------------------------------------------------------------
# Theophylline model using linCmt
m1 <- function() {
  ini({
    tka <- .5
    tcl <- -3.2
    tv <- -1
    eta.ka ~ 1
    eta.cl ~ 2
    eta.v ~ 1
    add.err <- 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.err)
  })
}

fit1 <- nlmixr(m1, theo_sd, est="saem", table=tableControl(cwres=TRUE, npde=TRUE))

print(fit1)

#fit1 <- fit1 %>% addCwres() # In case this was not specified under model fit prior to estimation one could add this here to the results

#GoF by xpose
xpdb <- xpose_data_nlmixr(fit1)

dv_vs_pred(xpdb) +
  ylab("Observed Theophylline Concentrations (ng/mL)") +
  xlab("Population Predicted Theophylline Concentrations (ng/mL)")
dv_vs_ipred(xpdb) +
  ylab("Observed Theophylline Concentrations (ug/mL)") +
  xlab("Individual Predicted Theophylline Concentrations (ng/mL)")
res_vs_pred(xpdb) +
  ylab("Conditional Weighted Residuals") +
  xlab("Population Predicted Theophylline Concentrations (ng/mL)")
res_vs_idv(xpdb) +
  ylab("Conditional Weighted Residuals") +
  xlab("Time (h)")
prm_vs_iteration(xpdb)
absval_res_vs_idv(xpdb, res = 'IWRES') +
  ylab("Individual Weighted Residuals") +
  xlab("Time (h)")
absval_res_vs_pred(xpdb, res = 'IWRES')  +
  ylab("Individual Weighted Residuals") +
  xlab("Population Predicted Theophylline Concentrations (ng/mL)")
ind_plots(xpdb, nrow=3, ncol=4) +
  ylab("Predicted and Observed Theophylline concentrations (ng/mL)") +
  xlab("Time (h)")
res_distrib(xpdb) +
  ylab("Density") +
  xlab("Conditional Weighted Residuals")
nlmixr::vpc(fit1,nsim=500, show=list(obs_dv=T),
            ylab = "Theophylline Concentrations (ng/mL)", xlab = "Time (h)")
