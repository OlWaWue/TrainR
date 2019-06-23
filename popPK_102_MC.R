##
## Funktion dient nur dazu einen Kreis zu malen
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}


## Funktion zum simulieren der Zufallswürfe
mc.pi = function(n, draw=T) {
  
  ## runif ist ein Zufallswert zwischen -1 und 1  
  ## wobei jeder Wert gleich Wahrscheinlich ist! NICHT wie rnorm => Werte um Mittelwert wahrscheinlicher!
  x = runif(n, -1, 1)
  y = runif(n, -1, 1)
  
  ## Berechnen ob die mit x und y bezeichnete Punkte im Kreis liegen => 1 oder nicht => 0
  pin = (ifelse(sqrt(x^2 + y^2 <= 1), 1, 0))
  
  ## Erzeuge datenframe fürs plotten
  more_dat <- data.frame(x, y, pin)
  
  ## Plotte wenn gewünscht
  if (draw==TRUE)
  {
    dat <- circleFun(c(0,0),2,npoints = 100)
    
    pl <- ggplot(dat,aes(x,y)) + geom_path() + geom_point(data=more_dat, aes(x=x,y=y, colour=factor(pin) )) +
      xlim(-1,1)+ylim(-1,1) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
            legend.position = "none") +
      geom_hline(aes(yintercept=-1)) +geom_hline(aes(yintercept=1)) +
      geom_vline(aes(xintercept=-1)) + geom_vline(aes(xintercept=1))
    print(pl)
    
  }
  ## Gibt pi angenähert zurück => 4-fache Warscheinlichkeit dass der Punkt 
  ## im Kreis liegt im Verhältnis zu allen Punkten
  return(4 * sum(pin)/n)
}

set.seed(1111)

mc.pi(10, T)



