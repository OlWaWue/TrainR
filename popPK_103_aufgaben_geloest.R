library(ggplot2)

V_pop <- 50 #L
omega_vc <- 0.05 # SD des Random Effects auf V

ke_pop <- 0.045 #1/h
omega_ke <- 0.12 # SD des Random Effects auf ke

ka_pop <- 1.25 # 1/h
omega_ka <- 1.34 # SD des Random Effects auf ka

F_oral = 0.8
AMT = 200
sim_times=seq(from=0, to=12, by=0.1)
tau=12

pk_model_ss<- function(AMT, F_oral, Vc, ke, ka, time, tau=12){
  conc <- F_oral*AMT/Vc *(ka/(ka-ke))*( (exp(-ke*(time) )/(1-exp(-ke*tau))) - (exp(-ka*(time) )/(1-exp(-ka*tau))) )
  return(conc)
}

obj_func_no_re <- function(tdm_data, par) {
  measured_conc = tdm_data$measured_conc
  simulated_conc = pk_model_ss(AMT, F_oral, par[1], par[2], par[3], tdm_data$time, tau)
  return(sum( (measured_conc-simulated_conc)^2 ))
}


tdm_data = data.frame(time=c(5,10), measured_conc=c(7.25, 4.75))

simulated_data = data.frame(time=sim_times, 
                            simulated_conc=pk_model_ss(AMT, F_oral, V_pop, ke_pop, ka_pop, sim_times, tau))

ggplot(simulated_data) + geom_line(aes(x=time, y=simulated_conc)) + 
  geom_point(data=tdm_data, aes(x=time, y=measured_conc))                            


optim_res <- optim(par=c(V_pop, ke_pop, ka_pop), fn=obj_func_no_re, tdm_data=tdm_data)


simulated_ind_data = data.frame(time=sim_times, 
                                simulated_conc=pk_model_ss(AMT, F_oral, optim_res$par[1], optim_res$par[2], optim_res$par[3], sim_times, tau))


ggplot(simulated_data) + geom_line(aes(x=time, y=simulated_conc)) + 
  geom_point(data=tdm_data, aes(x=time, y=measured_conc)) +
  geom_line(data=simulated_ind_data, aes(x=time, y=simulated_conc), colour="red")
