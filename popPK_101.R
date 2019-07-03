
library(ggplot2)
library(readxl) #diese library benötigt man zum Laden der Excel Datei

## Daten aus Excel-Datei laden

examplocin <- read_excel("examplocin.xlsx")

## Daten anschauen; head zeigt nur die ersten Zeilen an; tail würde die letzen Zeilen anzeigen

head(examplocin)

### Daten als Punkte in einer Abbildung darstellen; DV= dependent variable (z.B. Plasmakonzentration)

ggplot(examplocin) + geom_point(aes(x=TIME, y=DV))

### Punkte nach Patienten ID einfärben

ggplot(examplocin) + geom_point(aes(x=TIME, y=DV, colour=as.factor(ID) )) 

### Abbildungen nach geschlecht unterteilen

ggplot(examplocin) + geom_point(aes(x=TIME, y=DV, colour=as.factor(ID) )) + facet_wrap(facets = "SEX")

## und umgekehrt

ggplot(examplocin) + geom_point(aes(x=TIME, y=DV, colour=as.factor(SEX) )) + facet_wrap(facets = "ID") +
  geom_line(aes(x=TIME, y=DV))

## und mit halblogarithmischer Darstellung

ggplot(examplocin) + geom_point(aes(x=TIME, y=DV, colour=as.factor(SEX) )) + facet_wrap(facets = "ID") +
  geom_line(aes(x=TIME, y=DV)) + scale_y_log10()


## PK-Modell erstellen

pk_1cmt_bolus <- function(time, amt, Vd, Cl){
  ke <- Cl/Vd
  
  y <- amt/Vd*exp(-ke*time) 
  
  return(y)
}

## PK - Modell ausprobieren

times <- 0:24

sim_cp <- pk_1cmt_bolus(times, 100, 80, 15)

sim_pk_data <- data.frame(TIME = times, DV = sim_cp)

plot(sim_pk_data)

## Populationsmodell berechnen, welche Werte für Vd und Cl passsen im Schnitt am besten zu den 4 Patienten-Datensätzen

## TIME wird er aus den Daten nehmen, Vd und Cl wird er nicht in den Examplocin-Daten finden und weiß daher, 
## dass er diese beiden Parameter mittels nicht linearer Regression ermitteln soll

pk_mod_pop <- nls(DV~pk_1cmt_bolus(TIME, 100, Vd, Cl), start=list(Vd=80, Cl=15), data=examplocin )

## Ergebniss des popPK Modells

summary(pk_mod_pop)

## Aus dem Populationsmodell heraus simulieren, TIME findet er wieder in den Daten,
## Vd und Cl nimmt er aus dem Modell

examplocin$DV_pop_pred <- predict(pk_mod_pop, examplocin)

## Darstellung der Populationsvorhersage

ggplot(examplocin) + geom_point(aes(x=TIME, y=DV, colour=as.factor(SEX) )) + facet_wrap(facets = "ID") +
  geom_line(aes(x=TIME, y=DV_pop_pred))

## Individuelle Modelle für jeden Patienten erstellen
## Beachten => der Datensatz wird auf den jeweiligen Patienten beschränkt!!

pk_mod_1 <- nls(DV~pk_1cmt_bolus(TIME, 100, Vd, Cl), 
                start=list(Vd=80, Cl=15), 
                data=examplocin[examplocin$ID == 1,]) ## <-- hier nur Pat 1

pk_mod_2 <- nls(DV~pk_1cmt_bolus(TIME, 100, Vd, Cl), 
                start=list(Vd=80, Cl=15), 
                data=examplocin[examplocin$ID == 2,]) ## <-- hier nur Pat 2

pk_mod_3 <- nls(DV~pk_1cmt_bolus(TIME, 100, Vd, Cl), 
                start=list(Vd=80, Cl=15), 
                data=examplocin[examplocin$ID == 3,]) ## <-- hier nur Pat 3

pk_mod_4 <- nls(DV~pk_1cmt_bolus(TIME, 100, Vd, Cl), 
                start=list(Vd=80, Cl=15), 
                data=examplocin[examplocin$ID == 4,]) ## <- hier nur Pat 4


## Für jeden Patienten seine individuelle Vorhersage mit seinem eigenen Modell erstellen

examplocin$DV_ind_pred <- 0 ## <-- neue Spalte im Datensatz einfügen und alle Werte auf 0 setzen

examplocin$DV_ind_pred <- ifelse(examplocin$ID==1, ## <-- nur für Datensätze deren ID = 1 ist
                                 predict(pk_mod_1, examplocin[examplocin$ID == 1,]),  ## <-- simuliere Wert mit Modell für Pat 1
                                 examplocin$DV_ind_pred) ## <-- sonst bleibt der Wert, so wie er vorher war


## und das ganze für jeden Patienten

examplocin$DV_ind_pred <- ifelse(examplocin$ID==2, 
                                 predict(pk_mod_2, examplocin[examplocin$ID == 2,]), 
                                 examplocin$DV_ind_pred)

examplocin$DV_ind_pred <- ifelse(examplocin$ID==3, 
                                 predict(pk_mod_3, examplocin[examplocin$ID == 3,]), 
                                 examplocin$DV_ind_pred)

examplocin$DV_ind_pred <- ifelse(examplocin$ID==4, 
                                 predict(pk_mod_4, examplocin[examplocin$ID == 4,]), 
                                 examplocin$DV_ind_pred)

## Jetzt fügen wir die individuelle Vorhersage noch in die Abbildung mit ein

ggplot(examplocin) + geom_point(aes(x=TIME, y=DV, colour=as.factor(SEX) )) + facet_wrap(facets = "ID") +
  geom_line(aes(x=TIME, y=DV_pop_pred)) + ## <-- das die Populationsvorhersage
  geom_line(aes(x=TIME, y=DV_ind_pred), colour="red") ## <-- das ist die individuelle Vorhersage


## Berechnen wir die Abweichung der individuellen Vorhersage vom gemessenen Wert

## Residuelle Variabilität

examplocin$EPSILON <- examplocin$DV - examplocin$DV_ind_pred

summary(examplocin$EPSILON)

sd(examplocin$EPSILON)

ggplot(examplocin) + geom_point(aes(x=DV, y=DV_ind_pred)) + ## <-- gemessen vs. individuelle Vorhersage
  geom_abline(slope = 1, intercept = 0) + ## <-- diagonale = Idealfall
  xlim(0,1.5) + ylim(0,1.5) ## Limits der Achsen manuell festlegen

ggplot(examplocin) + geom_point(aes(x=DV, y=DV_pop_pred)) + ## <-- gemessen vs. Populationsvorhersage
  geom_abline(slope = 1, intercept = 0) + ## <-- diagonale = Idealfall
  xlim(0,1.5) + ylim(0,1.5)

## Untersuchen der einzelnen PK Modelle, welche Parameter gelten für welchen Patienten?

summary(pk_mod_pop)
summary(pk_mod_1)
summary(pk_mod_2)
summary(pk_mod_3)
summary(pk_mod_4)

## Alle PK Parameter in einer Tabelle zusammenfassen


pk_params <- rbind(coef(pk_mod_pop)) ## <- rbind fügt eine weitere Zeile hinzu
pk_params <- rbind(pk_params, coef(pk_mod_1)) ## <-- "coef" zieht die Koeffizienten aus dem Modell
pk_params <- rbind(pk_params, coef(pk_mod_2))
pk_params <- rbind(pk_params, coef(pk_mod_3))
pk_params <- rbind(pk_params, coef(pk_mod_4))

pk_params <- cbind(pk_params, ID=c(0, 1,2,3,4)) ## <- cbind fügt eine Spalte hinzu

pk_params <- cbind(pk_params, ETA_VD=pk_params[,1]-pk_params[1,1]) ## <- weitere Spalte enthält Differenz von ind_Vd und pop_Vd

pk_params <- cbind(pk_params, ETA_CL=pk_params[,2]-pk_params[1,2]) ## <- weitere Spalte enthält Differenz von ind_Cl und pop_Cl

pk_params <- as.data.frame(pk_params)

pk_params


### Jetzt fügen wir zu den PK Parametern Patientendaten dazu
## Vielleicht wird es so möglich die Abweichungen zu erklären?


examplocin$WT

pk_params$WT <- c( ((52+63+124+79)/4), 52,63,124,79)

examplocin$SEX

pk_params$SEX <- c(-1, 0,1,0,1)

pk_params


###
## -- Die Visuelle Auswertung zeigt ...
###


ggplot(pk_params[-1,]) + geom_point(aes(x=WT, y=ETA_VD)) ## <-- Ein Zusammenhang zwischen Vd und WT scheint naheliegend
 
ggplot(pk_params[-1,]) + geom_point(aes(x=WT, y=ETA_CL)) ## <-- Cl und WT scheinen nicht sehr offensichtlich zu korrelieren

ggplot(pk_params[-1,]) + geom_boxplot(aes(group=SEX, y=ETA_CL)) ## Offenbar besteht ein Zusammenhang zwischen Cl und SEX

ggplot(pk_params[-1,]) + geom_point(aes(x=ETA_VD, y=ETA_CL)) ## Vd und Cl scheinen nicht korreliert


## Genauer bekommt man eine Antwort mit der linearen regression
## die funktion für lineares modell "lm" ist in der "summary" funktion "eingenistet"

## Ohne das Ergebnis der "lm" funktion zu speichern wird die Zusammenfassung erstellt


summary(lm(WT~ETA_VD, data = pk_params[-1,]))

summary(lm(WT~ETA_CL, data = pk_params[-1,]))

summary(lm(SEX~ETA_CL, data = pk_params[-1,]))

summary(lm(ETA_VD~ETA_CL, data = pk_params[-1,]))



###
## -- PK Modell mit Covariaten
###

pk_mod_pop

pk_params

pk_1cmt_bolus_cov <- function(time, amt, THETA1, THETA2, WT, SEX){
  
  Vd_pop <- THETA1
  Vd_ind <- Vd_pop *(WT/79.5)
  
  Cl_pop <- THETA2
  Cl_ind <- Cl_pop * (1.2^SEX)
  
  ke <- Cl_ind/Vd_ind
  
  y <- amt/Vd_ind*exp(-ke*time)
  
  return(y)
}

## umcodieren der geschlechtswerte in den Daten m

examplocin$SEX <- ifelse(examplocin$SEX=="m", 1, 0)


pk_mod_pop_cov <- nls(DV~pk_1cmt_bolus_cov(TIME, 100, THETA1, THETA2, WT, SEX), start=list(THETA1=1, THETA2=1), data=examplocin )

## Ergebniss des popPK Modells

summary(pk_mod_pop_cov)

## Aus dem Populationsmodell mit Kovariation heraus simulieren, TIME, WT und SEX findet er wieder in den Daten,
## THETA1 und THETA2 im Modell

examplocin$DV_pop_pred_cov <- predict(pk_mod_pop_cov, examplocin)

## Darstellung der Populationsvorhersage mit Kovariatenmodell

ggplot(examplocin) + geom_point(aes(x=TIME, y=DV )) + facet_wrap(facets = "ID") +
  geom_line(aes(x=TIME, y=DV_pop_pred)) + geom_line(aes(x=TIME, y=DV_pop_pred_cov), colour="red")
