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
                     "SD for ETA on Vc:",
                     min = 0,
                     max = 2,
                     value = .1, step=0.01),
        numericInput("ke",
                     "elimination rate [1/h]:",
                     max=0.45, value=0.45, step=0.01),
        sliderInput("ke_SD",
                    "SD for ETA on ke:",
                    min = 0,
                    max = 2,
                    value = 0.1, step=0.01),
        numericInput("ka",
                     "absorption rate [1/h]:",
                     min=0.5, value=1.15, max=1.5, step = 0.1),
        sliderInput("ka_SD",
                    "SD for ETA on ka:",
                    min = 0,
                    max = 2,
                    value = 0.1, step=0.01),
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
                      "Steady state?", value=F),
        numericInput("tau",
                     "Interdose Interval [h]",
                     value=12, min=2, max=24),
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
    pk_parameter = NULL
  )
  
   output$mainPlot <- renderPlot({
     
     times <- seq(from=0,
                  to=input$tau,
                  by=0.1)
     

       
       ## Zufällig für jeden PK-Parameter input$n_mc Proben ziehen  
         
      mc_Vc <- input$Vc * exp(rnorm(n = input$n_mc, mean=0, sd=input$Vc_SD))
      mc_ke <- input$ke * exp(rnorm(n = input$n_mc, mean=0, sd=input$ke_SD))
      mc_ka <- input$ka * exp(rnorm(n = input$n_mc, mean=0, sd=input$ka_SD))
         ##
         ## Ablegen der simulierten Parameter in app_data
      app_data$pk_parameter <- data.frame(n=1:input$n_mc, Vc=mc_Vc, ke=mc_ke, ka=mc_ka)
         

         
      
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
                                     times,
                                     input$tau) 
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
       geom_vline(aes(xintercept=(plot_dat$TIME[plot_dat$DELTA == max(plot_dat$DELTA)])), colour="red") +
       annotate(geom = "label", x=(plot_dat$TIME[plot_dat$DELTA == max(plot_dat$DELTA)]), y=max(plot_dat$DELTA)/2, label=paste((plot_dat$TIME[plot_dat$DELTA == max(plot_dat$DELTA)]), "h"), size=8)+
       xlab("Time [h]") + ylab("Delta plasmaconcentration [mg/L]")
   })
   
   
   ## Zeichne die simulierten häufigkeitsverteilungen in einem dritten Plot
   output$histPlot <- renderPlot(({
     
     ## Nur wenn Daten im app_data$parameter vorhanden sind, tue etwas
     if(is.null(app_data$pk_parameter)){
       return()
     } else{
       #
       #Erzeuge je PK-Parameter einen Plot => Histogramm und Dichtekurve der zufällig erzeugten PK-Parameter
       pl1 <- ggplot(app_data$pk_parameter) +theme_bw() + geom_histogram(aes(x=Vc, y=..density..), bins=50, fill="white", colour="black") + 
         xlab("Value for Vc") + ylab("Density") +
         geom_density(aes(x=Vc, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)
       
       pl2 <- ggplot(app_data$pk_parameter) +theme_bw()+ geom_histogram(aes(x=ke, y=..density..), bins=50, fill="white", colour="black") + 
         xlab("Value for ke") + ylab("Density")+
         geom_density(aes(x=ke, y=..density..), colour="red", size=.5, fill="red",alpha=0.25, linetype=1)
       
       pl3 <- ggplot(app_data$pk_parameter) +theme_bw()+ geom_histogram(aes(x=ka, y=..density..), bins=50, fill="white", colour="black") + 
         xlab("Value for ka") + ylab("Density")+
         geom_density(aes(x=ka, y=..density..), colour="red", size=.5, fill="red", alpha=0.25, linetype=1)
       
       ## Stelle alle drei Plots nebeneinander in einer Zeile dar
       gridExtra::grid.arrange(pl1,pl2,pl3, nrow=1, ncol=3)
     }

   }))
}

# Run the application 
shinyApp(ui = ui, server = server)

