
library(ggplot2)

##
## Wichtig ist also die Variabilität in den PK - Parametern
## 

### Modelliert wird das allerdings nicht so wie wir das hier gesehen haben mit Mittelwert +/- SD, weil die 
### Parameter dabei theoretisch kleiner 0 werden können.

## Eine Variante, wie die Parameter stets über Null liegen ist die Darstellung zum Beispiel für Vc:

## Vc = Vpop * exp(ETA_VC)

## Wobei ETA_VC aus einer Normalverteilung mit dem Mittelwert 0 und der SD = omega kommt

## in R =>

Vpop <- 50 #L

omega_vc <- 0.05 # SD des Random Effects

ETA_VC <- rnorm(n=1000, mean=0, sd=omega_vc)

Vc <- Vpop * exp(ETA_VC)

dat <- data.frame(n=1:1000, ETA_VC, Vc)

pl1 <- ggplot(data = dat) + geom_histogram(aes(x=ETA_VC, y=..density..), fill="white", colour="black") + 
  geom_density(aes(y=..density.., x=ETA_VC)) + theme_bw()

pl1

pl2 <- ggplot(data = dat) + geom_histogram(aes(x=Vc, y=..density..), fill="white", colour="black") + 
  geom_density(aes(y=..density.., x=Vc)) + theme_bw()

pl2

##
## Da man nicht nur einen PK Parameter mit Unsicherheit hat, sondern immer mehrere können wir für die anderen Parameter auch eine
## ähnliche Formel aufstellen

ke_pop <- 0.045 #1/h

omega_ke <- 0.12

ETA_ke <- rnorm(n=1000, mean=0, sd=omega_ke)

ke <- ke_pop * exp(ETA_ke)

dat <- cbind(dat, ke)


## und für ka das gleiche

ka_pop <- 1.25 # 1/h

omega_ka <- 1.34

ETA_ka <- rnorm(n=1000, mean = 0, sd=omega_ka)

ka <- ka_pop * exp(ETA_ka)

dat <- cbind(dat, ka)



## Jetzt schauen wir uns die weiteren plots an

pl3 <- ggplot(data = dat) + geom_histogram(aes(x=ETA_ke, y=..density..), fill="white", colour="black") + 
  geom_density(aes(y=..density.., x=ETA_ke)) + theme_bw()



pl4 <- ggplot(data = dat) + geom_histogram(aes(x=ke, y=..density..), fill="white", colour="black") + 
  geom_density(aes(y=..density.., x=ke)) + theme_bw()


pl5 <- ggplot(data = dat) + geom_histogram(aes(x=ETA_ka, y=..density..), fill="white", colour="black") + 
  geom_density(aes(y=..density.., x=ETA_ka)) + theme_bw()


pl6 <- ggplot(data = dat) + geom_histogram(aes(x=ka, y=..density..), fill="white", colour="black") + 
  geom_density(aes(y=..density.., x=ka)) + theme_bw()


gridExtra::grid.arrange(pl1, pl2, pl3, pl4, pl5, pl6, nrow=3, ncol=2)

##
## Um die zufälligen Effekte alle zusammen zu halten bildet man die sogenannte omega-Matrix eine symmetrische Matrix deren 
## Diagonalen Elemente normalerweise die Varianz der Random Effects tragen
## in R =>

omega_mat <- matrix(c(omega_vc^2, 0, 0,
                      0, omega_ke^2, 0,
                      0, 0, omega_ka^2), nrow = 3, ncol=3)


omega_mat

# In die Nicht-Diagnolen Felder stehen die Kovarianzen (wenn zwei random Effekte miteinander korrelieren)
#
##
## Lasst uns schauen, ob die zufälligen Effekte korrelieren

pl7 <- ggplot(data = dat) + geom_point(aes(x=ETA_ka, y=ETA_ke)) + theme_bw()
pl8 <- ggplot(data = dat) + geom_point(aes(x=ETA_ka, y=ETA_VC)) + theme_bw()
pl9 <- ggplot(data = dat) + geom_point(aes(x=ETA_ke, y=ETA_VC)) + theme_bw()

gridExtra::grid.arrange(pl7, pl3+coord_flip(), pl5+scale_y_reverse(),  nrow=2, ncol=2)
gridExtra::grid.arrange(pl8, pl1+coord_flip(), pl5+scale_y_reverse(),  nrow=2, ncol=2)
gridExtra::grid.arrange(pl9, pl1+coord_flip(), pl3+scale_y_reverse(),  nrow=2, ncol=2)

##
## -- puhh --- zum Glück keine Korrelationen zu erkennen. Die "off-diagnol" Felder bleiben bei Null ...... erstmal :-)
##


## Erzeugen wir unsere Simulationen


## PK model for steady state 1-cmt oral first order absorption and first order elimination
pk_model_ss<- function(AMT, F_oral, Vc, ke, ka, time, tau=12){
  conc <- F_oral*AMT/Vc *(ka/(ka-ke))*( (exp(-ke*(time) )/(1-exp(-ke*tau))) - (exp(-ka*(time) )/(1-exp(-ka*tau))) )
  return(conc)
}

times = seq(from=0, to=24, by=0.1)

## 1000 Kurven mit diesem Modell generieren

dat_mc <- NULL

  ## MC simulation der PK-Resultate n Verschiedene Kurven erzeugen je nachdem ob steady state oder single dose
  for(i in 1:1000){

      CP_mc <-  pk_model_ss(AMT=200, 
                            F_oral=0.87, 
                            dat$Vc[i], 
                            dat$ke[i], 
                            dat$ka[i], 
                            times) 
    ## Jede Spalte enthält eine Simulation => 1000 Spalten werden im laufe der Schleife erzeugt
    dat_mc <- cbind(dat_mc, CP_mc)
  }

## Bewirkt ein drehen der Tabelle um 90° => Spalten werden zu Zeilen und Zeilen zu spalten
## In jeder Spalte stehen jetzt die 1000 Simulationen für jeden einzelnen Zeitpunkt
## Jede Spalte steht für einen Zeitpunkt der simuliert wurde (diese Zeitpunkte sind in times hinterlegt)
transposed_data <- t(dat_mc)

## Jede Spalte wird jetzt auf 11 Zeilen eingeschrumpft => Eine Zeile für das 5, 10, 15, 20, 25% Percentil, Median 
## und 75,80,85,90 und 95% Percentil
s <- apply(transposed_data,2,function(x) quantile(x,probs=c(0.05,0.1,0.15,0.2,0.25,0.5,0.75,0.8,0.85,0.9,0.95)))

# diese Perzentile liegen danach in s[1, ], s[2, ],  s[3, ] ....

## Ergebnisse zusammenfassen in einem data.frame
plot_dat <- data.frame(TIME=times, CP_min_1=s[1,], CP_min_2=s[2,], CP_min_3=s[3,], CP_min_4=s[4,], CP_min_5=s[5,], 
                       CP=s[6,], 
                       CP_max_1=s[7,], CP_max_2=s[8,],CP_max_3=s[9,],CP_max_4=s[10,],CP_max_5=s[11,])

my_beautiful_colour = "red"

pl <- ggplot(plot_dat) + geom_line(aes(x=TIME, y=CP), colour=my_beautiful_colour, size=0.75) + 
  geom_ribbon(aes(x=TIME, ymax=CP_max_5, ymin=CP_min_1), alpha=0.15, fill=my_beautiful_colour) +
  geom_ribbon(aes(x=TIME, ymax=CP_max_4, ymin=CP_min_2), alpha=0.15, fill=my_beautiful_colour) +
  geom_ribbon(aes(x=TIME, ymax=CP_max_3, ymin=CP_min_3), alpha=0.15, fill=my_beautiful_colour) +
  geom_ribbon(aes(x=TIME, ymax=CP_max_2, ymin=CP_min_4), alpha=0.15, fill=my_beautiful_colour) +
  geom_ribbon(aes(x=TIME, ymax=CP_max_1, ymin=CP_min_5), alpha=0.15, fill=my_beautiful_colour) +
  theme_bw() + xlab("Time [h]") + ylab("Plasmaconcentration [mg/L]")

pl


