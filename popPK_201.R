library(mrgsolve)
library(tidyr)
library(dplyr)
library(ggplot2)

library(nlmixr)
library(xpose)
library(xpose.nlmixr)

mod <- mread("Examplotinib")

set.seed(16111983)

patient_data <- data.frame(ID=1:20, SEX=rbinom(20, 1, 0.5), WT=rnorm(20, 75, 15), CLCR=rnorm(20,115, 35))



res <- mod %>% ev(ID=1:20, amt=1000, addl=0, ii=24) %>% idata_set(patient_data) %>% carry_out(amt,evid, addl, ii) %>% mrgsim(end=24)

plot(res)

dat <- as.data.frame(res)



sample_dat <- filter(dat, dat$time %in% c(0, 0.5, 2, 4, 8, 12, 24))

del_col <- which(colnames(sample_dat) %in% c("GUT", "CENT", "PER", "KA_ind", "CL_ind", "V_ind", "KE_ind", "VP_ind", "Q_ind", "IPRED"))

data <- sample_dat[,-del_col]

data$cmt <- 2

data$mdv <- 0

data$cmt <- ifelse(data$evid == 1, 1, 2)

data$mdv <- ifelse(data$evid == 1, 1, 0)

data$LNDV <- log(data$DV)




data <- data[-which(data$time==0 & data$evid==0),]


write.csv(data, "pat_dat.csv", row.names = F)

ggplot(data, aes(x=time, y=DV)) + geom_point() + facet_wrap(data$ID)



ggplot(data, aes(x=time, y=DV, color=as.factor(ID) )) +
  geom_point() +
  geom_line()

ggplot(data, aes(x=time, y=DV)) +
  geom_point() +
  stat_smooth() +
  scale_y_log10()

my2CptModel <- function() {
  ini({
    tka <- log(0.22)
    tcl <- log(3.395)
    tv2  <- log(35)
    tv3  <- log(450)
    tq   <- log(5)
    eta.ka ~ 1
    eta.cl ~ 1
    eta.v2 ~ 1
    eta.v3 ~ 1
    eta.q ~ 1
    add.err <- 0.075
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v2 <- exp(tv2 + eta.v2)
    v3 <- exp(tv3 + eta.v3)
    q  <- exp(tq + eta.q)
    d/dt(depot) = -ka * depot
    d/dt(center) = ka * depot - cl / v2 * center + q/v3 * periph - q/v2 * center
    d/dt(periph) = q/v2 * center - q/v3 * periph
    cp = center / v2
    cp ~ add(add.err)
  })
}


my2CptFit <- nlmixr(my2CptModel, data, est="saem", 
                    control=saemControl(nBurn = 200, nEm=300),
                    table=tableControl(cwres=TRUE, npde=TRUE))

plot(my2CptFit)


xpdb <- xpose.nlmixr::xpose_data_nlmixr(my2CptFit)




dv_vs_pred(xpdb) +
  ylab("Observed Examplotinib Concentrations (mg/L)") +
  xlab("Population Predicted Examplotinib Concentrations (mg/L)")
dv_vs_ipred(xpdb) +
  ylab("Observed Examplotinib Concentrations (mg/L)") +
  xlab("Individual Predicted Examplotinib Concentrations (mg/L)")
res_vs_pred(xpdb) +
  ylab("Conditional Weighted Residuals") +
  xlab("Population Predicted Examplotinib Concentrations (mg/L)")
res_vs_idv(xpdb) +
  ylab("Conditional Weighted Residuals") +
  xlab("Time (h)")
prm_vs_iteration(xpdb)
absval_res_vs_idv(xpdb, res = 'IWRES') +
  ylab("Individual Weighted Residuals") +
  xlab("Time (h)")
absval_res_vs_pred(xpdb, res = 'IWRES')  +
  ylab("Individual Weighted Residuals") +
  xlab("Population Predicted Examplotinib Concentrations (mg/L)")
ind_plots(xpdb, nrow=3, ncol=4) +
  ylab("Predicted and Observed Examplotinib concentrations (mg/L)") +
  xlab("Time (h)")
res_distrib(xpdb) +
  ylab("Density") +
  xlab("Conditional Weighted Residuals")
nlmixr::vpc(my2CptFit,nsim=500, show=list(obs_dv=T),
            ylab = "Examplotinib Concentrations (mg/L)", xlab = "Time (h)")
