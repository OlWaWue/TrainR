## Bibliothek für graphische Darstellung laden
library(ggplot2)

## Daten aus dem R allgemein Paket laden
data("Theoph")

## PK Daten für Theophyllin anschauen
Theoph

# Wt = Körpergewicht in kg
# Dose = mg/kg Körpergewicht
# Time in Stunden
# conc in mg/L


## Vorbereitung: Berechnen wir die Dosis in mg für jeden Patienten

Theoph$oral_dose <- Theoph$Dose*Theoph$______

##
## Aufgabe 1: Die Daten betrachten: Vervollständige den Code unten
##


# a) Alle Daten in einem Scatterplot (Punktewolke) darstellen Conc. vs. Time
ggplot(Theoph) + geom_point(aes(x=______, y=______))


# b) Punkte so einfärben, dass jedes Individuum eine eigene Farbe bekommt
ggplot(Theoph) + geom_point(aes(x=______, y=______, colour=______))

# c) Punkte in ihrer Größe Verändern, dass Punktgröße das Körpergewicht widerspiegelt
ggplot(Theoph) + geom_point(aes(x=______, y=______, colour=______, size=______))

# d) Die Punkte mit linien verbinden
ggplot(Theoph) + geom_point(aes(x=______, y=______, colour=______, size=______)) +
  ______(aes(x=______, y=______, colour=______))

# e) Jedes Individuum in einem eigenen Graphen darstellen

ggplot(Theoph) + geom_point(aes(x=______, y=______, colour=______, size=______)) +
  ______(aes(x=______, y=______, colour=______)) + _____(facets = Theoph$Subject)

# f) Daten in einer halblogarithmischen Darstellung anschauen

ggplot(Theoph) + geom_point(aes(x=______, y=______, colour=______, size=______)) +
  ______(aes(x=______, y=______, colour=______)) + _____(facets = Theoph$Subject) + _____()


##
## Prima! Welches PK Modell passt wahrscheinlich am besten? 
##

##
## Aufgabe 2: Schreibe eine Funktion mit Kinetik erster Ordnung absorption und Elimination erster Ordnung
## Parameter für die Funktion sind Dosis, F_oral, Vd, ke, ka und time
##

# Definiere die Funktion

pk_model <- function(Dosis, F_oral, Vd, ke, ka, time){
  conc <- _______
  return(conc)
}

# Erstelle eine Zeitsequenz (Stunden von 0 bis 24 in 0.1h Abstand)
time_sequence <- seq(from=0, to=24, by=0.1)

# Berechne für diese Zeiten die Konzentration für Dosis=500mg, F_oral=0.8, Vd=30L, ke=0.05 1/h, und ka=0.9 1/h

concentrations <- pk_model(_____)

# Verbinde die Daten zu einer Tabelle, einem data.frame

sim_pk_data <- data.frame(Time=time_sequence, conc=concentrations)

# Schau Dir die Tabelle an

# nur die ersten sechs Zeilen
head(_____)

# nur die letzten sechs Zeilen
tail(_____)

# die Struktur des data.frames (Welche Art Daten liegen darin?)
str(_____)


##
## Verwende diese Daten um daraus einen Graphen zu erstellen
##

# a) Stelle conc vs Time als Linie dar

ggplot(sim_pk_data) + _____(aes(x=Time, y=conc))


# b)  Füge eine Beschritung der x und y Achse ein

ggplot(sim_pk_data) + geom_line(aes(x=Time, y=conc)) + xlab("Zeit [h]") + ____("Plasmakonzentration [mg/L]")

# c) Nun mit halblogarithmischer Darstellung

ggplot(sim_pk_data) + geom_line(aes(x=Time, y=conc)) + xlab("Zeit [h]") + ylab("Plasmakonzentration [mg/L]") +____


## Versuche die populations-Parameter für Theophyllin
## Verwende den non-linear least squares Algorithmus um die PK-Parameter zu finden, die auf alle Patienten zutreffen

pk_parm_pop <- nls(conc~____(Dosis = oral_dose,
                              F_oral = 0.8,
                              Vd,
                              ke,
                              ka,
                              time = Time), start=list(Vd=30, ke=0.05, ka=0.9), data=Theoph )


## Schauen wir uns die popPK Parameter an:

summary(pk_parm_pop)

# Nun Versuchen wir die Daten vom ersten Patienten zu isolieren um für diesen Patienten die PK-Parameter zu berechnen

theoph_patient_1 <- Theoph[Theoph$Subject==____,]


# Berechnen die individuellen Parameter für Patient 1, Startbedingungen sollte poPK Parameter werden
pk_parm_1 <- nls(conc~____(Dosis = oral_dose,
                                 F_oral = 0.8,
                                 Vd,
                                 ke,
                                 ka,
                                 time = Time), start=list(Vd=___, ke=___, ka=___), data=theoph_patient_1 )

summary(pk_parm_1)

# Wir machen dies jetzt für jeden einzelnen Patienten mit Hilfe einer for-schleife

# Für 1 bis 12 Patienten
for(i in 1:max(Theoph$Subject)){
  # nimm aktuellen Patient vor
  theoph_patient_cur <- Theoph[Theoph$Subject==___,]
  
  
  # Berechnen die individuellen Parameter für Patient 1, Startbedingungen sollte poPK Parameter werden
  pk_parm_cur <- nls(conc~pk_model(Dosis = oral_dose,
                                 F_oral = 0.8,
                                 Vd,
                                 ke,
                                 ka,
                                 time = Time), start=list(Vd=26.3, ke=0.08, ka=1.55), data=_______ )
  
  print(summary(pk_parm_cur))
}

## Anstatt die Daten für jeden Auszudrucken sollen alle individuellen Parameter in einer Tabelle gespeichert werden
##
## Dafür definieren wir die Tabelle

pk_params = data.frame()


# Für 1 bis 12 Patienten
for(i in 1:max(Theoph$Subject)){
  
  _______
  
  ## Ziehe die aktuellen Parameter raus
  sum <- summary(pk_parm_cur)
  Vd_cur <- sum$coefficients[1,1]
  ke_cur <- sum$coefficients[1,2]
  ka_cur <- sum$coefficients[1,3]
  
  cur_data <- data.frame(ID=i, WT=Theoph[Theoph$Subject==i,]$Wt[1], Vd=Vd_cur, ke=ke_cur, ka=ka_cur)
  pk_params <- rbind(pk_params, cur_data)
}

pk_params


##
## Gibt es einen Zusammenhang eines PK-Parameters und Gewicht?

plot(pk_params$WT~pk_params$___)
abline(lm(pk_params$WT~pk_params$____))
