library(deSolve)

SIR_model<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+I+R)
    lam<-parameters['beta']*I/P
    
    dS <- mu*P-mu*S-lam*S         
    dI <- lam*S-tau*I-mu*I 
    dR <- tau*I-mu*R
    list(c(dS, dI, dR))}) 
}


initP<-999 
initI<-1 
initR<-0 
initS<-initP-initI-initR
istate <- c(S = initS, I = initI, R = initR)

parameters <- c(
  mu=(1/(50*52*7)),    
  beta=0.25,              
  tau=1/10            
)

R0<-parameters['beta']/
  (parameters['tau']+parameters['mu'])

time_start <- 0
time_stop <- 365 
deltat<- 0.1
tps <- seq(time_start , time_stop , by = deltat)

out <- deSolve::ode(y = istate, times = tps, func = SIR_model, parms = parameters) 

Flu2009<-as.data.frame(out)
# plot(Flu2009)
plot(Flu2009$S,Flu2009$I,ylim = c(0,1000), xlim = c(0,1000))


SIR_model_dd<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+I+R)
    lam<-parameters['beta']*I
    
    dS <- mu*P-mu*S-lam*S         
    dI <- lam*S-tau*I-mu*I 
    dR <- tau*I-mu*R
    list(c(dS, dI, dR))}) 
}
out_dd <- deSolve::ode(y = istate, times = tps, func = SIR_model_dd, parms = parameters) 

Flu2009_dd<-as.data.frame(out_dd)
lines(Flu2009_dd$S,Flu2009_dd$I, col="sandybrown")

plot(Flu2009_dd$I)

################# different parameter sets   #############################################
parameters_pp <- c(
  mu=(1/(50*52*7)),   #### change the values inside this vector
  beta=0.25,              
  tau=1/10            
)
parameters_p <- c(
  mu=(1/(50*52*7)),   #### change the values inside this vector 
  beta=0.25,              
  tau=1/10            
)
parameters_mm <- c(
  mu=(1/(50*52*7)),   #### change the values inside this vector 
  beta=0.25,              
  tau=1/10            
)
parameters_m <- c(
  mu=(1/(50*52*7)),  #### change the values inside this vector  
  beta=0.25,              
  tau=1/10            
)

#################   store different solutions ###############################################

### model solutions with different parameters
out_pp <- deSolve::ode(y = istate, times = tps, func = SIR_model_dd, parms = parameters_pp) 
out_p <- deSolve::ode(y = istate, times = tps, func = SIR_model_dd, parms = parameters_p) 

out_mm <- deSolve::ode(y = istate, times = tps, func = SIR_model_dd, parms = parameters_mm) 
out_m <- deSolve::ode(y = istate, times = tps, func = SIR_model_dd, parms = parameters_m) 

### model solutions with different models
out_interventions <- deSolve::ode(y = istate, times = tps, func = SIR_model_interventios, parms = parameters) 

## model solutions with different models and parameters
out_interventions_pp <- deSolve::ode(y = istate, times = tps, func = SIR_model_interventios, parms = parameters_pp) 


##############  adding lines to plots   ##########################################################

plot(out[,3], ylim=c(0,1000))
plot(out[,3])






