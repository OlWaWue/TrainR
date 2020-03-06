
data("Theoph")

Theoph

## Dosis ändern

Theoph$Dose <- Theoph$Dose*Theoph$Wt

library(ggplot2)

ggplot(data=Theoph, mapping=aes(x=Time, y=conc)) + geom_point()

## Einzelne Individuen rausfinden

ggplot(data=Theoph, mapping=aes(x=Time, y=conc, colour=Subject)) + geom_point()

## Alternative: facet_wrap

ggplot(data=Theoph, mapping=aes(x=Time, y=conc)) + geom_point() + facet_wrap(facets="Subject")


## Nur ein Individuum rauspicken

data_sub_1 <- Theoph[Theoph$Subject==1,]

print(data_sub_1)

ggplot(data=data_sub_1, mapping=aes(x=Time, y=conc)) + geom_point()

## Modell schreiben

pk_model <- function(AMT, TIME, KA, KE, VD, BIOV) {
  
  conc <- BIOV*AMT/VD * KA/(KA-KE) * (exp(-KE*TIME)-exp(-KA*TIME))
  
  return(conc)
}

data_sub_1$PRED <- pk_model(AMT=data_sub_1$Dose, TIME=data_sub_1$Time, KA=1, KE=0.5, VD=30, BIOV=1)

ggplot(data=data_sub_1, mapping=aes(x=Time, y=conc)) + geom_point() + geom_line(mapping=aes(y=PRED))

## orale Bioverfügbarkeit von Theophyllin ist 0.76%

objective_function <- function(par=c(KA=1.5, KE=0.5, VD=30), ACTUAL_DATA, AMT, BIOV){
  
  measured_conc <- ACTUAL_DATA$conc
  
  pred_conc <- pk_model(AMT=ACTUAL_DATA$Dose, 
                        TIME=ACTUAL_DATA$Time, 
                        KA=par[["KA"]], 
                        KE=par[["KE"]], 
                        VD=par[["VD"]], 
                        BIOV)
  

  return(sum( (measured_conc-pred_conc)^2 ))
}

## optimalen Wert finden

optim_res_ls <- optim(par=c(KA=1.5, KE=0.5, VD=30), fn=objective_function, 
                      ACTUAL_DATA=data_sub_1, BIOV=0.76)

optim_res_ls$par


data_sub_1$PRED <- pk_model(AMT=data_sub_1$Dose, TIME=data_sub_1$Time, 
                            KA=optim_res_ls$par[["KA"]], 
                            KE=optim_res_ls$par[["KE"]], 
                            VD=optim_res_ls$par[["VD"]], BIOV=0.76)

ggplot(data=data_sub_1, mapping=aes(x=Time, y=conc)) + geom_point() + geom_line(mapping=aes(y=PRED))

tail(Theoph)

for(i in 1:12){
  data_cur_sub <- Theoph[Theoph$Subject==i,]
  
  optim_res_curr <- optim(par=c(KA=1.5, KE=0.5, VD=30), fn=objective_function, 
                          ACTUAL_DATA=data_cur_sub, BIOV=0.76)
  
  print(optim_res_curr$par)
  
  
  data_cur_sub$PRED <- pk_model(AMT=data_cur_sub$Dose, TIME=data_cur_sub$Time, 
                              KA=optim_res_curr$par[["KA"]], 
                              KE=optim_res_curr$par[["KE"]], 
                              VD=optim_res_curr$par[["VD"]], BIOV=0.76)
  
  pl <- ggplot(data=data_cur_sub, mapping=aes(x=Time, y=conc)) + geom_point() + geom_line(mapping=aes(y=PRED))
  
  print(pl)
}

results <- NULL

for(i in 1:12){
  data_cur_sub <- Theoph[Theoph$Subject==i,]
  
  optim_res_curr <- optim(par=c(KA=1.5, KE=0.5, VD=30), fn=objective_function, 
                          ACTUAL_DATA=data_cur_sub, BIOV=0.76)
  
  
  new_ind_result <- data.frame(Subject=i, 
                               WT=data_cur_sub$Wt[1], 
                               KA=optim_res_curr$par[["KA"]], 
                               KE=optim_res_curr$par[["KE"]], 
                               VD=optim_res_curr$par[["VD"]])
  
  results <- rbind(results, new_ind_result) ##
  
  data_cur_sub$PRED <- pk_model(AMT=data_cur_sub$Dose, TIME=data_cur_sub$Time, 
                                KA=optim_res_curr$par[["KA"]], 
                                KE=optim_res_curr$par[["KE"]], 
                                VD=optim_res_curr$par[["VD"]], BIOV=0.76)
  
  pl <- ggplot(data=data_cur_sub, mapping=aes(x=Time, y=conc)) + geom_point() + geom_line(mapping=aes(y=PRED))
  
  print(pl)
}

summary(results)

ggplot(results, aes(x=WT, y=VD)) + geom_point() + geom_smooth(method="lm")


objective_function_stat_mod <- function(par=c(ETA_KA=0, ETA_KE=0, ETA_VD=0), KA, KE, VD, ACTUAL_DATA, AMT, BIOV){
  
  measured_conc <- ACTUAL_DATA$conc
  
  pred_conc <- pk_model(AMT=ACTUAL_DATA$Dose, 
                        TIME=ACTUAL_DATA$Time, 
                        KA=KA+par[["ETA_KA"]], 
                        KE=KE+par[["ETA_KE"]], 
                        VD=VD+par[["ETA_VD"]], 
                        BIOV)
  
  
  return(sum( (measured_conc-pred_conc)^2 ))
}

results <- NULL

for(i in 1:12){
  data_cur_sub <- Theoph[Theoph$Subject==i,]
  
  optim_res_curr <- optim(par=c(ETA_KA=0, ETA_KE=0, ETA_VD=0), fn=objective_function_stat_mod, 
                          ACTUAL_DATA=data_cur_sub, BIOV=0.76, KA=2.19, KE=0.089, VD=24.26)
  
  
  new_ind_result <- data.frame(Subject=i, 
                               WT=data_cur_sub$Wt[1], 
                               ETA_KA=optim_res_curr$par[["ETA_KA"]], 
                               ETA_KE=optim_res_curr$par[["ETA_KE"]], 
                               ETA_VD=optim_res_curr$par[["ETA_VD"]])
  
  results <- rbind(results, new_ind_result) ##
  
  data_cur_sub$PRED <- pk_model(AMT=data_cur_sub$Dose, TIME=data_cur_sub$Time, 
                                KA=2.19+optim_res_curr$par[["ETA_KA"]], 
                                KE=0.089+optim_res_curr$par[["ETA_KE"]], 
                                VD=24.26+optim_res_curr$par[["ETA_VD"]], BIOV=0.76)
  
  pl <- ggplot(data=data_cur_sub, mapping=aes(x=Time, y=conc)) + geom_point() + geom_line(mapping=aes(y=PRED))
  
  print(pl)
}

summary(results)


