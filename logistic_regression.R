library(ggplot2)
library(tidyverse)
library(dplyr)
library(pROC)
library(broom)

my_dat <- data.frame(HbA1c=runif(100, 5,9), Neuropathy=rbinom(100, 1, 0.5))

my_dat$HbA1c <- ifelse(my_dat$Neuropathy==1, my_dat$HbA1c*1.2, my_dat$HbA1c)


ggplot(my_dat) + geom_point(aes(y=Neuropathy, x=HbA1c))

my_mod <- glm(data = my_dat, formula = Neuropathy~HbA1c, family=binomial)

new_dat <- data.frame(HbA1c = seq(3,15,1))

new_dat$pred <- predict(my_mod, newdata = new_dat, type="response")


ggplot(my_dat) + geom_point(aes(y=Neuropathy, x=HbA1c)) + geom_line(data= new_dat, aes(y=pred, x=HbA1c))

glance(my_mod)

my_roc <- roc(response=my_dat$Neuropathy, predictor = my_mod$fitted.values)

my_roc

plot(my_roc)

coords(my_roc, "best")
