library(survival)

mydata <- lung

mydata$recodedStatus <- ifelse(mydata$status==1, 0, 1)

head(mydata)

mysurv <- Surv(time=mydata$time, event=mydata$recodedStatus)

myfit <-survfit(mysurv~mydata$sex)

survdiff(mysurv~mydata$sex)

myfit


plot(myfit, col = c("blue", "red"), conf.int = T, mark=3)

