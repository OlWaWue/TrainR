##
## -- Shiny application for Teaching 
##

library(shiny)
library(ggplot2)

## PK model for single dose 1-cmt oral first order absorption and first order elimination
pk_model <- function(Dosis, F_oral, Vd, ke, ka, time){
  conc <- F_oral*Dosis/Vd *(ka/(ka-ke))*(exp(-ke*time)-exp(-ka*time))
  return(conc)
}


## PK model for steady state 1-cmt oral first order absorption and first order elimination
pk_model_ss<- function(Dosis, F_oral, Vd, ke, ka, time, tau=12){
  conc <- F_oral*Dosis/Vd *(ka/(ka-ke))*( (exp(-ke*(time) )/(1-exp(-ke*tau))) - (exp(-ka*(time) )/(1-exp(-ka*tau))) )
  return(conc)
}

# Define UserInterface for application 
ui <- fluidPage(
   
   # Application title
   titlePanel("SimPK - by Oliver Scherf-Clavel - JMU Wuerzburg 2019"),
   
   # Sidebar containing input Elements
   sidebarLayout(
      sidebarPanel(
        numericInput("Vc",
                     "central Volume [L]:",
                     min=1, value=30),
         sliderInput("Vc_SD",
                     "SD for Vc [L]:",
                     min = 0,
                     max = 5,
                     value = 0, step=0.1),
        numericInput("ke",
                     "elimination rate [1/h]:",
                     max=0.45, value=0.45, step=0.01),
        sliderInput("ke_SD",
                    "SD for ke [1/h]:",
                    min = 0,
                    max = 0.1,
                    value = 0, step=0.01),
        numericInput("ka",
                     "absorption rate [1/h]:",
                     min=0.5, value=1.15, max=1.5, step = 0.1),
        sliderInput("ka_SD",
                    "SD for ka [1/h]:",
                    min = 0,
                    max = 0.2,
                    value = 0, step=0.01),
        numericInput("AMT",
                     "Dose [mg]:",
                     min=1, value=5000),
        numericInput("F_oral",
                     "systemically available fraction:",
                     min=0, max=1,
                     value=0.8),
        checkboxInput("logTrans",
                      "y-axis in log-10", value=F),
        checkboxInput("isSS",
                      "Steady state? tau=12 h", value=F),
        checkboxInput("do_mc",
                      "Monte Carlo Simulation?", value=F),
        numericInput("n_mc",
                     "Number of MC Samples",
                     min=1, max=10000, value=1000, step=100)

      ),
      
      # Show a plot of the generated PK-Profile
      mainPanel(
         plotOutput("mainPlot"),
         plotOutput("secPlot"),
         plotOutput("histPlot")
      )
   )
)

# Define server logic - Computations behind the surface
server <- function(input, output) {
   
  app_data <- reactiveValues(
    plot_data = NULL,
    pk_parameter = NULL,
    pk_par_dens = NULL
  )
  
   output$mainPlot <- renderPlot({
     
     times <- seq(from=0,
                  to=12,
                  by=0.1)
     
     ## Je nachdem ob angehakt ist nimm das steady state oder das single dose modell
     if(input$isSS){ 
       
       ## Berechne Kurve für eingetragenen Werte
       CP <- pk_model_ss(input$AMT, 
                      input$F_oral, 
                      input$Vc, 
                      input$ke, 
                      input$ka, times)
       
       ## Berechne Kurve für PK-Parameter 2 SD nach unten abweichen
       CP_min <- pk_model_ss(input$AMT, 
                          input$F_oral, 
                          (input$Vc-2*input$Vc_SD), 
                          (input$ke-2*input$ke_SD), 
                          (input$ka-2*input$ka_SD), 
                          times)
       
       ## Berechne Kurve für PK-Parameter 2 SD nach oben abweichen
       CP_max <- pk_model_ss(input$AMT, 
                          input$F_oral, 
                          (input$Vc+2*input$Vc_SD), 
                          (input$ke+2*input$ke_SD), 
                          (input$ka+2*input$ka_SD), 
                          times) 
       
       
     } else {
         CP <- pk_model(input$AMT, 
                    input$F_oral, 
                    input$Vc, 
                    input$ke, 
                    input$ka, times)
     
          CP_min <- pk_model(input$AMT, 
                    input$F_oral, 
                    (input$Vc-2*input$Vc_SD), 
                    (input$ke-2*input$ke_SD), 
                    (input$ka-2*input$ka_SD), 
                    times)
          CP_max <- pk_model(input$AMT, 
                        input$F_oral, 
                        (input$Vc+2*input$Vc_SD), 
                        (input$ke+2*input$ke_SD), 
                        (input$ka+2*input$ka_SD), 
                        times) 
       }
     
     plot_dat <- data.frame(TIME=times, CP=CP, CP_min=CP_min, CP_max=CP_max, DELTA=(CP_max-CP_min))
     
     ## Was ist hier nicht abzubilden? => Der Fall wenn ein Wert nach oben und einer nach unten abweicht! 
     ## => Monte-Carlo-Simulation
     ## => man nimmt z.B. 1000 verschiedene Kurven auf, wobei man zufällig 1000 Kombinationen aus Normalverteilungen zieht
     if (input$do_mc){ 
       
       ## Zufällig für jeden PK-Parameter input$n_mc Proben ziehen aus 
       ## der Verteilung mit Modalwert: input$Vc/ke/ka und Standardabweichung input$VC_SD/ke_SD/ka_SD
         
         mc_Vc <- rnorm(n = input$n_mc, mean=input$Vc, sd=input$Vc_SD)
         mc_ke <- rnorm(n = input$n_mc, mean=input$ke, sd=input$ke_SD)
         mc_ka <- rnorm(n = input$n_mc, mean=input$ka, sd=input$ka_SD)
         ##
         ## Ablegen der simulierten Parameter in app_data
         app_data$pk_parameter <- data.frame(n=1:input$n_mc, Vc=mc_Vc, ke=mc_ke, ka=mc_ka)
         
         ##
         ## Das Folgende dient nur der Darstellung der Gaußkurve über den gezogenen Proben 
         ## => Abgleich ob die gezogenen Proben normalverteilt sind
         ##
         
         ##
         ## Es wird die Dichtefunktion der Normalverteilung benutzt um über den Bereich vom kleinssten zum größten Parameter
         ## in jeder Gruppe die Dichtefunktion mit den gegebenen Mittelwert und SD zu erzeugen
         
         dens_Vc <- data.frame(x=seq(from=min(mc_Vc), to=max(mc_Vc), by=0.001) , 
                                     y= dnorm(x=seq(from=min(mc_Vc), to=max(mc_Vc), by=0.001), mean=input$Vc, sd=input$Vc_SD ) )
         
         dens_ke <- data.frame(x=seq(from=min(mc_ke), to=max(mc_ke), by=0.001) , 
                                     y= dnorm(x=seq(from=min(mc_ke), to=max(mc_ke), by=0.001), mean=input$ke, sd=input$ke_SD ) )
       
         dens_ka <- data.frame(x=seq(from=min(mc_ka), to=max(mc_ka), by=0.001) , 
                                     y= dnorm(x=seq(from=min(mc_ka), to=max(mc_ka), by=0.001), mean=input$ka, sd=input$ka_SD ) )
         
         ##
         ## Ablegen der Dichtekurven in app_data als Liste. Die pk_par_dens ist also eine Liste aus data.frames also eine Liste aus
         ## Tabellen, kann man sich als Schubladen vorstellen. Schublade 1, 2 und 3 und in jeder Schublade liegt eine
         ## data.frame Tabelle
         app_data$pk_par_dens <- list(dens_Vc, dens_ke, dens_ka)
         
      
         ## temporäre Variable zum Sammeln der Simulationen erzeugen
         dat_mc <- NULL
         
         ## Fortschrittsanzeige benutzen
         withProgress(message = "Performing Monte Carlo Simulation", max = input$n_mc, {
         
           ## MC simulation der PK-Resultate n Verschiedene Kurven erzeugen je nachdem ob steady state oder single dose
             for(i in 1:input$n_mc){
               if (input$isSS){
                  CP_mc <-  pk_model_ss(input$AMT, 
                                     input$F_oral, 
                                     mc_Vc[i], 
                                     mc_ke[i], 
                                     mc_ka[i], 
                                     times) 
               } else {
                 CP_mc <-  pk_model(input$AMT, 
                                       input$F_oral, 
                                       mc_Vc[i], 
                                       mc_ke[i], 
                                       mc_ka[i], 
                                       times) 
               }
               dat_mc <- cbind(dat_mc, CP_mc)
               incProgress(1)
             }
         })
         
        s <- apply(t(dat_mc),2,function(x) quantile(x,probs=c(0.1,0.5,0.9)))
        
        ## Überschreibe plot_dat mit den Monte Carlo Resultaten
        plot_dat <- data.frame(TIME=times,CP_min=s[1,],CP=s[2,],CP_max=s[3,], DELTA=(s[3,]-s[1,]))

     } else {
       ##
       ## Frühere Simulationsergebnisse löschen, wenn das Häkchen entfernt wird
       app_data$pk_parameter <- NULL
       app_data$pk_par_dens <- NULL
     }
     
     

     ## Ablegen der Plot Daten in app_data
     app_data$plot_data <- plot_dat
     
     ## Erzeuge PK-Kurve
     pl <- ggplot(plot_dat) + geom_line(aes(x=TIME, y=CP)) + 
       geom_ribbon(aes(x=TIME, ymax=CP_max, ymin=CP_min), alpha=0.15, fill="red") +
       geom_line(aes(x=TIME, y=CP_min), colour="red") + 
       geom_line(aes(x=TIME, y=CP_max), colour="red") + theme_bw() + xlab("Time [h]") + ylab("Plasmaconcentration [mg/L]")
     
     ## Wenn Häkchen bei logTrans => Element der Log-Transformation hinzufügen
     if(input$logTrans)
       pl <- pl + scale_y_log10()
     
     print(pl)
     
   })
   
   output$secPlot <- renderPlot({
     #
     # Lade Plot-Daten aus app_data
     plot_dat <- app_data$plot_data
     
     ## Erzeuge den Sensitivitätsplot
     ggplot(plot_dat) + geom_line(aes(x=TIME, y=DELTA)) + theme_bw() + 
       xlab("Time [h]") + ylab("Delta plasmaconcentration [mg/L]")
   })
   
   
   ## Zeichne die simulierten häufigkeitsverteilungen in einem dritten Plot
   output$histPlot <- renderPlot(({
     
     ## Nur wenn Daten im app_data$parameter vorhanden sind, tue etwas
     if(is.null(app_data$pk_parameter)){
       return()
     } else{
       #
       #Erzeuge je PK-Parameter einen Plot => Histogramm die zufällig erzeugten PK-Parameter
       # => Kurve ist die per Dichtefunktion erzeugte ideale Verlaufsform für jeden PK-Parameter. Daten ist eine Liste von data.frames 
       ## Elemente in Listen kann man mit [[Eintrag Nummer]] aufrufen. An Stelle 1 steht der data.frame für Vc, an 2 für ke und an 3 für ka
       ## Weil in dieser Reihenfolge die Elemente in die Liste eingefügt wurden (siehe oben)
       pl1 <- ggplot(app_data$pk_parameter) +theme_bw() + geom_histogram(aes(x=Vc, y=..density..), fill="white", colour="black") + 
         xlab("Value for Vc") + ylab("Density") +
         geom_line(data=app_data$pk_par_dens[[1]], aes(x=x, y=y), colour="blue", size=1.25, alpha=0.5, linetype=1)
       
       pl2 <- ggplot(app_data$pk_parameter) +theme_bw()+ geom_histogram(aes(x=ke, y=..density..), fill="white", colour="black") + 
         xlab("Value for ke") + ylab("Density")+
         geom_line(data=app_data$pk_par_dens[[2]], aes(x=x, y=y), colour="blue", size=1.25, alpha=0.5, linetype=1)
       
       pl3 <- ggplot(app_data$pk_parameter) +theme_bw()+ geom_histogram(aes(x=ka, y=..density..), fill="white", colour="black") + 
         xlab("Value for ka") + ylab("Density")+
         geom_line(data=app_data$pk_par_dens[[3]], aes(x=x, y=y), colour="blue", size=1.25, alpha=0.5, linetype=1)
       
       ## Stelle alle drei Plots nebeneinander in einer Zeile dar
       gridExtra::grid.arrange(pl1,pl2,pl3, nrow=1, ncol=3)
     }

   }))
}

# Run the application 
shinyApp(ui = ui, server = server)

