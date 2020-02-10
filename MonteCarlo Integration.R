library(ggplot2)



function_to_integrate <- function(x,pars = c(ka = 2,
                                             ke=0.55,
                                             Vd=50,
                                             D=500,
                                             f_abs = 0.6)) {
  ka <- pars[['ka']]
  ke <- pars[['ke']]
  Vd <- pars[['Vd']]
  D <- pars[['D']]
  f_abs <- pars[['f_abs']]
  y <- D*f_abs/Vd * ka/(ka-ke) * (exp(-ke*x) - exp(-ka*x))
  return(y)
}


xlab = "Time since dose"
xunit = "h"

ylab = "Plasmaconcentration"
yunit = "mg/L"

## Integrate between
xmin=0
xmax=24


function_data <- data.frame(x=seq(xmin, xmax, 0.1))
function_data$y <- function_to_integrate(function_data$x)


ymin=0
ymax=max(function_data$y)

n_samples <- 10000

x_total <- runif(n=n_samples, min=xmin, max=xmax)
y_total <- runif(n=n_samples, min=ymin, max=ymax)


x_function <- ifelse(y_total <= function_to_integrate(x_total), x_total, NA)
y_function <- ifelse(is.na(x_function), NA, y_total)
  

mc_int_data <- data.frame(x_total, y_total, x_function, y_function)

mc_int_data$is_under_curve <- ifelse(is.na(x_function), F, T)

points_under_the_curve <- sum(!is.na(x_function))

total_area <- (xmax-xmin)*(ymax-ymin)

area_under_the_curve <- points_under_the_curve*total_area/n_samples

ggplot(mc_int_data) + geom_point(aes(x=x_total, y=y_total, colour=is_under_curve)) +
  ggtitle(paste("Area under the Curve:", round(area_under_the_curve,2), xunit, "x", yunit ),
          paste("Total points: ", n_samples, " - ", "Points under the curve: ", points_under_the_curve, " (", round(points_under_the_curve/n_samples,2)*100, " %)", sep="" )) + 
  theme_bw() +
  geom_line(data=function_data, aes(x=x, y=y), size=0.8)+
  xlab(paste(xlab, "[", xunit, "]")) + ylab(paste(ylab, "[", yunit, "]")) +
    scale_colour_manual(values=c("TRUE"="blue",
                                 "FALSE"="firebrick"),
                        guide = guide_legend(override.aes = list(
                          linetype =  c(0,0),
                          fill = c("blue", "firebrick")
                        ),
                        title="Is under the curve?"), 
                        labels=c("FALSE", 
                                 "TRUE")) 
