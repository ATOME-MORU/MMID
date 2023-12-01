library(deSolve)
setwd('C:/mgh - mmid')

SEIR_model<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+E+I+R)
    P_aux<-S+E+max(rnorm(1,I,30),0)+R
    lam<- max(rnorm(1,parameters['beta'],5e-3),0)*max(rnorm(1,I,30),0)/P_aux
    
    dS <- mu*P-mu*S-lam*S   
    dE <- lam*S-gamma*E-mu*E
    dI <- gamma*E-tau*I-mu*I 
    dR <- tau*I-mu*R
    
    list(c(dS, dE, dI, dR))}) 
}


initP<-1000
initE<-0
initI<-1 
initR<-0 
initS<-initP-initI-initR-initE
istate <- c(S = initS, E = initE, I = initI, R = initR)

parameters <- c(
  mu=(1/(50*52*7)),    
  beta=.25,              
  gamma=1/4,          
  tau=1/10            
)

R0<-parameters["beta"]*parameters["gamma"]/
  ((parameters["gamma"]+parameters["mu"])*(parameters["tau"]+parameters["mu"]))
R0
time_start <- 0
time_stop <- 365 
deltat<- 0.01
tps <- seq(time_start , time_stop , by = deltat)

out <- deSolve::ode(y = istate, times = tps, func = SEIR_model, parms = parameters) 

Flu2009<-as.data.frame(out)
  
plot(out)
# plot(Flu2009)

N<-rowSums(Flu2009[,2:5])
plot(N)

library(EpiEstim)
Flu2009$incidence<-parameters["gamma"]*Flu2009$E
ww<-which(Flu2009$time %% 1 ==0.0)

res_parametric_si <- estimate_R(abs(Flu2009$incidence[ww]), method="parametric_si", config = make_config(list( mean_si = 9.0, std_si = 1.5)) ) 
plot(res_parametric_si, legend = FALSE)

res_parametric_si$R



SEIR_model_4n<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+E+I1+I2+I3+I4+R)
    n<-4
    lam<- rnorm(1,parameters['beta'],5e-3)*(I1+I2+I3+I4)/P
    dS <- mu*P-mu*S-lam*S   
    dE <- lam*S-gamma*E-mu*E
    dI1 <- gamma*E-tau*n*I1-mu*I1 
    dI2 <- tau*n*I1-mu*I2-tau*n*I2
    dI3 <- tau*n*I2-mu*I3-tau*n*I3
    dI4 <- tau*n*I3-mu*I4-tau*n*I4
    dR <- n*tau*I4-mu*R
    list(c(dS, dE, dI1, dI2, dI3, dI4, dR))}) 
}

initP<-1000
initE<-0
initI1<-1 
initI2<-1 
initI3<-1 
initI4<-1 
initR<-0 
initS<-initP-initI1-initI2-initI3-initI4-initR-initE
istate4n <- c(S = initS, E = initE, I1 = initI1, I2 = initI2, I3 = initI3, I4 = initI4, R = initR)

out4n <- deSolve::ode(y = istate4n, times = tps, func = SEIR_model_4n, parms = parameters) 

Flu2009_4n<-as.data.frame(out4n)
plot(out4n)
plot(Flu2009_4n)

