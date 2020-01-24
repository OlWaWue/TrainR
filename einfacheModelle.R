
library(ggplot2)

## Zum Beispiel eine Kalibriergerade
## Messungen per UV-Absorption (y)
## gegen die Konzentration (x)
data_for_lin_reg <- data.frame(Absorption=c(0.04, 0.12, 0.19, 0.30, 0.35, 0.48),
                               Konzentration=c(0, 0.1, 0.2, 0.3, 0.4, 0.5))


print(data_for_lin_reg)


ggplot(data = data_for_lin_reg, aes(x=Konzentration, y=Absorption)) + geom_point() +
  geom_smooth(method = "lm") # Fügt regressionsgerade hinzu

## Unser Modell besteht hier aus der Annahme, dass die Absorption proportional mit der Konzentraiton
## zunimmt. Das Modell besteht also aus der Steigung und die durch Eigenabsorption
## des Lösungsmittels verursachten y-Achsenabschnitt. 
## Die einzige unabhängige Variable hier ist die Konzentration, die 
## Modellparameter m (Steigung) und b (Achsenabschnitt) werden im Vektor 'par' 
## zusammengefasst

## => f(Konzentration|par) = m * Konzentration + b

## Das mathematische Modell wird in R mit Hilfe einer 'function' definiert

model <- function(Konzentration, par) {
  
  m <- par[1]
  b <- par[2]
  
  return(m*Konzentration+b) ## Vorhersage der Absorption anhand der Konzentration
  
}

## Wie finden wir die Modellparameter für die 'beste' Geradengleichung heraus?
## => Optimieren der Objektiven Funktion
## Die Objektive Funktion stellt ein Maß dafür dar, wie 'gut' die gemessenen Daten zu einem
## parameterset par (THETA) passen

obj_func <- function(Daten, par) {
  
  Absorption <- Daten$Absorption         ## y1,...,i
  Konzentration <- Daten$Konzentration   ## x1,...,i
  
  ob_func_value <- sum( (Absorption-model(Konzentration,par))^2 ) # Quadrierte Abweichung von gemessener Absorption 
                                                                  # soll minimiert werden
                                                                  # damit negative differenzen auch 'bestraft' werden
                                                                  # => Quadrieren
  return(ob_func_value)
}     


## Jetzt rufen wir R auf um die obj_func zu minimieren
## => Finde den Satz von parametern 'par', der die Quadrierte Abweichung von Gemessenem zu Berechnetem Wert 
## minimiert => Maximum-Likelihood-Schätzer für 'par'

fit <- optim(par = c(m=1,b=0),
             fn = obj_func,
             Daten=data_for_lin_reg)

par_hat <- fit$par


print(par_hat)

## mit diesen Parametern können wir nun die Regressionsgerade erstellen


data_with_lin_reg <- data_for_lin_reg

## Neue Spalte anfügen

data_with_lin_reg$berechnete_Absorption <- model(data_with_lin_reg$Konzentration, par_hat)


print(data_with_lin_reg)

ggplot(data = data_with_lin_reg, aes(x=Konzentration, y=Absorption)) + geom_point() +
  geom_smooth(method = "lm") + # Fügt automatische Regressionsgerade hinzu
  geom_line(aes(x=Konzentration, y=berechnete_Absorption), colour="red") 

## liegt genau auf der von R erstellten Regressionsgeraden!

## Wie genau bekommen wir allerdings die Konfidenzintervalle?
## Man benötigt die Standardfehler der geschätzten parameter 'par_hat'



## Dazu lassen wir uns praktischerweise die hessian Matrix gleich mitberechnen


fit_with_hessian <- optim(par = c(m=1,b=0),
                    fn = obj_func,
                    Daten=data_for_lin_reg,
                    hessian = TRUE) ## Berechne diese Matrix

print(fit_with_hessian$hessian)


## Wir benötigen außerdem die residuelle Varianz


## Berechnung der residuellen Varianz mittels Maximum likelihood (Zur Minimierung wurde OLS verwendet)
sig2 <- fit$value/(nrow(data_for_lin_reg)-2) ## Abzug von zwei (Modellparametern) ## von den Freiheitsgeraden

print(sig2)

var_cov_mat <- solve(fit_with_hessian$hessian/(2*sig2))

RSE <- sqrt(diag(var_cov_mat))/par_hat*100

RSE

SE <- sqrt(diag(var_cov_mat))

SE

CI_95_oberes <- par_hat + 1.96*sqrt(diag(var_cov_mat))
CI_95_unteres <- par_hat - 1.96*sqrt(diag(var_cov_mat))



print(paste("Ermittelte Steigung:", 
            round(par_hat[1], 5), 
            "(95% KI:", 
            round(CI_95_unteres[1], 5),
            "-", 
            round(CI_95_oberes[1], 5),
            "), SE:", 
            round(SE[1], 5)
            
            )
      )

print(paste("Ermittelter Achsenabschnitt:", 
            round(par_hat[2], 5), 
            "(95% KI:", 
            round(CI_95_unteres[2], 5),
            "-", 
            round(CI_95_oberes[2], 5),
            "), SE:", 
            round(SE[2], 5)
            
            )
      )

print(paste("Residueller Standardfehler:", round(sqrt(sig2),5) ))


## Vergleiche mit der bereits in R hinterlegten Methode zur
## linearen Regression

automatic_fit <- lm(data = data_for_lin_reg, formula = Absorption~Konzentration)

summary(automatic_fit)


## Hier die Konfidenzintervalle für die Grundgesamtheit

data_with_lin_reg_with_CI <- data_with_lin_reg

se_a <- function(x, x_2, mean, n, sig2) {
  
  return(sqrt(sig2*(1/n+(x-mean)^2/(x_2) )))
  
}



data_with_lin_reg_with_CI$CI_Absorption_low <- data_with_lin_reg_with_CI$berechnete_Absorption+se_a(data_with_lin_reg_with_CI$Konzentration, 
                                                                                                      mean(data_with_lin_reg_with_CI$Konzentration), 
                                                                                                      sum((data_with_lin_reg_with_CI$Konzentration-mean(data_with_lin_reg_with_CI$Konzentration))^2),
                                                                                                      nrow(data_with_lin_reg_with_CI), sig2) * qt(0.975, df=4)
data_with_lin_reg_with_CI$CI_Absorption_high <- data_with_lin_reg_with_CI$berechnete_Absorption- se_a(data_with_lin_reg_with_CI$Konzentration, 
                                                                                                       mean(data_with_lin_reg_with_CI$Konzentration), 
                                                                                                       sum((data_with_lin_reg_with_CI$Konzentration-mean(data_with_lin_reg_with_CI$Konzentration))^2),
                                                                                                       nrow(data_with_lin_reg_with_CI), sig2)* qt(0.975, df=4)


print(data_with_lin_reg_with_CI)

## Darstellung
ggplot(data = data_with_lin_reg_with_CI, aes(x=Konzentration, y=Absorption)) + geom_point() +
  geom_smooth(method = "lm") + # Fügt automatische Regressionsgerade hinzu
  geom_line(aes(x=Konzentration, y=berechnete_Absorption), colour="red") +
  geom_ribbon(aes(x=Konzentration, ymin=CI_Absorption_low, ymax=CI_Absorption_high), fill="red", alpha=0.25)



### Konfidenzintervall für zukünftige Beobachtungen


se_b <- function(x, x_2, mean, n, sig2) {
  
  return(sqrt(sig2*(1+1/n+(x-mean)^2/(x_2) )))
  
}



data_with_lin_reg_with_CI$CI_Absorption_low <- data_with_lin_reg_with_CI$berechnete_Absorption+se_b(data_with_lin_reg_with_CI$Konzentration, 
                                                                                                     mean(data_with_lin_reg_with_CI$Konzentration), 
                                                                                                     sum((data_with_lin_reg_with_CI$Konzentration-mean(data_with_lin_reg_with_CI$Konzentration))^2),
                                                                                                     nrow(data_with_lin_reg_with_CI), sig2) * qt(0.975, df=4)
data_with_lin_reg_with_CI$CI_Absorption_high <- data_with_lin_reg_with_CI$berechnete_Absorption- se_b(data_with_lin_reg_with_CI$Konzentration, 
                                                                                                       mean(data_with_lin_reg_with_CI$Konzentration), 
                                                                                                       sum((data_with_lin_reg_with_CI$Konzentration-mean(data_with_lin_reg_with_CI$Konzentration))^2),
                                                                                                       nrow(data_with_lin_reg_with_CI), sig2)* qt(0.975, df=4)


print(data_with_lin_reg_with_CI)

## Darstellung
ggplot(data = data_with_lin_reg_with_CI, aes(x=Konzentration, y=Absorption)) + geom_point() +
  geom_smooth(method = "lm") + # Fügt automatische Regressionsgerade hinzu
  geom_line(aes(x=Konzentration, y=berechnete_Absorption), colour="red") +
  geom_ribbon(aes(x=Konzentration, ymin=CI_Absorption_low, ymax=CI_Absorption_high), fill="red", alpha=0.25)



      