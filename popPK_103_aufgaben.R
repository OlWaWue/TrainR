library(ggplot2)

## Definition der PK-Parameter mit Variabilität

V_pop <- 50 #L
omega_vc <- 0.8 # SD des Random Effects auf V

ke_pop <- 0.045 #1/h
omega_ke <- 0.55 # SD des Random Effects auf ke

ka_pop <- 1.25 # 1/h
omega_ka <- 1.25 # SD des Random Effects auf ka

### Fix in dieser Übung
F_oral = 0.8
AMT = 200
sim_times=seq(from=0, to=12, by=0.1)
tau=12

pk_model_ss<- function(AMT, F_oral, Vc, ke, ka, time, tau=12){
  conc <- F_oral*AMT/Vc *(ka/(ka-ke))*( (exp(-ke*(time) )/(1-exp(-ke*tau))) - (exp(-ka*(time) )/(1-exp(-ka*tau))) )
  return(conc)
}

## Was macht das hier?
set.seed(11121987)

### bis hier fix



### par ist ein Vektor von Parametern c(V, ke, ka) der optimiert wird 
### damit das Ergebnis der OBJ-Fun möglichst klein wird
obj_func_ls <- function(tdm_data, par) {
  measured_conc = tdm_data$measured_conc
  simulated_conc = pk_model_ss(AMT, F_oral, par[1], par[2], par[3], tdm_data$time, tau)
  return(sum( (measured_conc-simulated_conc)^2 ))
}

### Diese Messdaten werden aufgenommen und in einer Tabelle abgelegt
tdm_data = data.frame(time=c(4, ## h post dose
                             7, ## h post dose
                             10), ## h post dose
                      measured_conc=c(9.25, #mg/L
                                      6.7,  #mg/L
                                      5.75)) #mg/L

head(tdm_data)

## Simuliere ein interdose interval mit popPK Parameter
simulated_data = data.frame(time=sim_times, 
                            simulated_conc=pk_model_ss(AMT, F_oral, V_pop, ke_pop, ka_pop, sim_times, tau))

## Abbildung erstellen um simulierte Konzentrationen als Linie und gemessene Konzentrationen als Punkte darzustellen
ggplot(simulated_data) + geom_line(aes(x=time, y=________)) + 
  geom_point(data=tdm_data, aes(x=time, y=________), size=3, shape=1) + theme_bw() +
  xlab("Zeit seit letzter Dosis [h]") + ylab("Konzentration [mg/L]")



### Jetzt versuchen wir die Kurve an die Messwerte anzupassen.
## Zunächst mit Least Squares 

## ls steht für least squares => Die Summe der Abweichungsquadrate ist minimal
## Die Funktion optim verändert die in "par" eingetragenen Parameter so lange, bis das Ergebnis der
## obj_func_ls minimal (also Abweichungsquadrate) minimal sind
optim_res_ls <- optim(par=c(V_pop, ke_pop, ka_pop), fn=obj_func_ls, tdm_data=tdm_data)

## Das sind die Werte für Vc, ke und ka, welche der Algorithmus errechnet hat
## NUR UNTER DER VORAUSSETZUNG, DASS ABWEICHUNGSQUADRATE MINIMAL SEIN SOLLEN
## Algorithmus hat also "freie Hand" was die Auswahl der PK-Parameter angeht
optim_res_ls$par

# Daten aufbereiten simulieren und gleich in ein data.frame packen
simulated_ind_data_ls = data.frame(time=sim_times, 
                                   simulated_conc=pk_model_ss(AMT, ## Dosis in mg
                                                           F_oral, ## oral systemisch verfügbarer Anteil
                                                           optim_res_ls$par[1], ## Vc
                                                           optim_res_ls$par[2], ## ke
                                                           optim_res_ls$par[3], ## ka
                                                           sim_times,  ## Zeiten zu simulieren
                                                           tau )) ## Interdose interval

### Erweitern um die "gefittete" Kurve
ggplot(simulated_data) + geom_line(aes(x=time, y=simulated_conc)) + 
  geom_point(data=tdm_data, aes(x=time, y=measured_conc), size=3, shape=1) +
  geom_line(data=________, aes(x=time, y=simulated_conc), colour="red") + theme_bw() +
  xlab("Zeit seit letzter Dosis [h]") + ylab("Konzentration [mg/L]")


#### Als nächstes wollen wir vorwissen aus dem PK-Modell nutzen!
## Der Algorithmus wird in der Wahl der PK-Parameter nicht mehr "alleine gelassen",
## wir geben ihm Grenzen vor, in denen er sich bewegen darf.

### par ist jetzt ein Vektor von ETAS c(ETA_V, ETA_ke, ETA_ka) der optimiert wird 
### damit das Ergebnis der OBJ-Fun möglichst klein wird
### !!! ES WERDEN JETZT ALSO DIE WAHRSCHEINLICHSTEN ETAS FÜR DIESEN PATIETEN GESUCHT !!!
### omegas ist jetzt ein Vektor, der die SD der random Effekte (der ETAs) enthält c(omega_Vc, omega_ke, omega_ka)
### pop-values ist ein Vektor der die populationsparameter beinhaltet
obj_func_bay <- function(tdm_data, par, omegas, pop_vals) {
  measured_conc = tdm_data$measured_conc
  Vc <- pop_vals[1] * exp(par[1]) ## Vc= V_pop x e^ETA_vc
  ke <- pop_vals[2] * exp(par[2]) ## ke= ke_pop x e^ETA_ke
  ka <- pop_vals[3] * exp(par[3]) ## ka= ka_pop x e^ETA_ka
  simulated_conc = pk_model_ss(AMT, F_oral,Vc, ke, ka, tdm_data$time, tau)
  
  ## Die Quadrate der Differenz sind noch da + Die gewichteten Quadrate der Abweichung ETAS von 0
  ## Es wird durch die Varianz gewichtet! 
  obj_res <- sum( (measured_conc-simulated_conc)^2 ) + sum( ((par-0)^2)/(omegas^2)) 
  return(obj_res)
}

### par beinhaltet jetzt die Startwerte für ETAs (0 ist das wahrscheinlichste ETA => modus der Verteilung)
optim_res_bay <- optim(par=c(0, 0, 0), 
                   fn=obj_func_bay, 
                   tdm_data=tdm_data, 
                   omegas=c(omega_vc, omega_ke, omega_ka),
                   pop_vals=c(V_pop, ke_pop, ka_pop))

### Die so ermitteln "optimalen" Parameter sind jetzt die, unter Einbeziehen des Vorwissens,
### für diesen Patienten WAHRSCHEINLICHSTEN ETAs
ind_etas <- optim_res_bay$par


## Aus den ETAs lassen sich nun die individuellen PK-Parameter berechnen
v_ind <- V_pop*exp(ind_etas[1]) # L
ke_ind <- ke_pop*exp(ind_etas[2]) #1/h
ka_ind <- ka_pop*exp(ind_etas[3]) # 1/h

## ls steht für least squares => Die Summe der Abweichungsquadrate ist minimal
simulated_ind_data_bay = data.frame(time=sim_times, 
                                   simulated_conc=pk_model_ss(AMT, F_oral, 
                                                              v_ind, 
                                                              ke_ind, 
                                                              ka_ind, sim_times, tau))

## Jetzt erweitern wir den Plot noch um die "bayesisch" gefitte Kurve
ggplot(simulated_data) + geom_line(aes(x=time, y=simulated_conc)) + 
  geom_point(data=tdm_data, aes(x=time, y=measured_conc), size=3, shape=1) +
  geom_line(data=simulated_ind_data_ls, aes(x=time, y=simulated_conc), colour="red", linetype=2) +
  geom_line(data=________, aes(x=time, y=simulated_conc), colour="red") + theme_bw() +
  xlab("Zeit seit letzter Dosis [h]") + ylab("Konzentration [mg/L]")



## Zum Abschluss möchten wir noch wissen, wo in der ETA-Verteilung der ganzen Population
## sich die individuellen ETAs diesen Patienten befinden
## Dazu simulieren wir 5000 ETAs für V ke und ka mit der jeweiligen SD omega_xx
MC_ETA_V <- rnorm(n=5000, mean=0, sd=________)
MC_ETA_ke <- rnorm(n=5000, mean=0, sd=________)
MC_ETA_ka <- rnorm(n=5000, mean=0, sd=________)

eta_dat <- data.frame(ETA_VC=MC_ETA_V,
                      ETA_KE=MC_ETA_ke,
                      ETA_KA=MC_ETA_ka )

## Abbildungen vorbereiten
## Als senkrechte (vertikale) Linie wird in jede Häufigkeitsverteilung das individuell ermittelte ETA eingepasst
pl_vc <- ggplot(data = eta_dat) + geom_histogram(aes(x=ETA_VC, y=..density..), fill="white", colour="black", bins=100) + 
         geom_density(aes(y=..density.., x=ETA_VC)) + theme_bw() +
         geom_vline(aes(xintercept=ind_etas[1]), size=3, alpha=0.5, colour="red") +
         xlab("Wert für ETA_Vc") + ylab("Häufigkeit")+ xlim(c(-5,5))

pl_ke <- ggplot(data = eta_dat) + geom_histogram(aes(x=ETA_KE, y=..density..), fill="white", colour="black", bins=100) + 
         geom_density(aes(y=..density.., x=ETA_KE)) +
         geom_vline(aes(xintercept=ind_etas[2]), size=3, alpha=0.5, colour="red") + theme_bw() +
         xlab("Wert für ETA_ke") + ylab("Häufigkeit")+ xlim(c(-5,5))

pl_ka <- ggplot(data = eta_dat) + geom_histogram(aes(x=ETA_KA, y=..density..), fill="white", colour="black", bins=100) + 
         geom_density(aes(y=..density.., x=ETA_KA)) +
         geom_vline(aes(xintercept=ind_etas[3]), size=3, alpha=0.5, colour="red") + theme_bw() +
         xlab("Wert für ETA_ka") + ylab("Häufigkeit") + xlim(c(-5,5))

## Alle Plots schon säuberlich untereinander ausgeben
gridExtra::grid.arrange(pl_vc, pl_ke, pl_ka, ncol=1, nrow=3)

##
## Aufgaben: Wie verändert sich die Simulation wenn die Standardabweichungen immer größer gesetzt werden 
## => Größenordnung 10 - 100 => sog. uninformed prior!
## 1. Was glaubt ihr?
## 2. was seht ihr tatsächlich?

## Wie verändert sich die Simulation wenn ein weiterer Messwert aufgenommen wurde time=0.5, c=8.5

