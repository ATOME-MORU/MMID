library(deSolve)
setwd('C:/mgh - mmid')

SACRD<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+A+C+R)
    lam<-parameters['beta']*(A+C)/P
    
    dS <- mu*P-mu*S-lam*S   
    dA <- (1-gamma)*lam*S - taua*A - mu*A
    dC <- gamma*lam*S - tauc*C - mu*C 
    dR <- (1-phi)*tauc*C + taua*A - mu*R
    dD <- phi*tauc*C
    list(c(dS, dA, dC, dR, dD))}) 
}


initP<-999 
initA<-0
initC<-1
initR<-0 
initD<-0
initS<-initP-initA-initC-initR
istate <- c(S = initS, A = initA, C = initC, R = initR, D = initD)

parameters <- c(
  mu = (1/(50*52*7)),    
  beta = 0.25, 
  gamma = 0.2,
  taua = 1/10,
  tauc = 1/10,
  phi=0.05
)

time_start <- 0
time_stop <- 600
deltat<- 0.1
tps <- seq(time_start , time_stop , by = deltat)

out <- deSolve::ode(y = istate, times = tps, func = SACRD, parms = parameters) 
Result<-as.data.frame(out)

plot(out)
plot(Result)
plot(Result$S,Result$C+Result$A, ylim = c(0,500),xlim = c(0,1000))


##########################################################################################################################

#                                   Compare with SIR_model

########################################################################################################################

SIR_model2<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+I+R)
    lam<-parameters['beta']*I/P
    
    dS <- mu*P-mu*S-lam*S         
    dI <- lam*S-(1-phi)*tau*I-delta*phi*I-mu*I 
    dR <- (1-phi)*tau*I-mu*R
    dD <- delta*phi*I
    list(c(dS, dI, dR, dD))}) 
}

initP<-999 
initI<-1 
initR<-0 
initD<-0
initS<-initP-initI-initR
istate2 <- c(S = initS, I = initI, R = initR, D = initD)

parameters2 <- c(
  mu = (1/(50*52*7)),    
  beta = 0.25,              
  tau = 1/10,
  delta = 1/10,
  phi=0.05
)

out2 <- deSolve::ode(y = istate2, times = tps, func = SIR_model2, parms = parameters2) 
Result2<-as.data.frame(out2)

lines(Result2$S,Result2$I, col="dodgerblue3",lwd=3)


#############    What if asymptomatic infections last longer?  ###################

parameters['taua']<-1/20
out3 <- deSolve::ode(y = istate, times = tps, func = SACRD, parms = parameters) 
Result3<-as.data.frame(out3)

lines(Result3$S,Result3$C+Result3$A, col="tomato3",lwd=3)
plot(Result3$C+Result3$A)
lines(Result3$S)


parameters['beta']<-0.09
parameters['taua']<-1/10
out_lowbeta <- deSolve::ode(y = istate, times = tps, func = SACRD, parms = parameters) 
parameters['taua']<-1/20
out3_lowbeta <- deSolve::ode(y = istate, times = tps, func = SACRD, parms = parameters) 

plot(out_lowbeta[,2],out_lowbeta[,3]+out_lowbeta[,4], ylim = c(0,150),xlim = c(0,1000))
lines(out3_lowbeta[,2],out3_lowbeta[,3]+out3_lowbeta[,4], col="seagreen2",lwd=3)


###########    Symptoms dictate infectiousness     ############################

SACRD_rho<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+A+C+R)
    lam<-parameters['beta']*(A*rho+C)/P
    
    dS <- mu*P-mu*S-lam*S   
    dA <- (1-gamma)*lam*S - taua*A - mu*A
    dC <- gamma*lam*S - tauc*C - mu*C 
    dR <- (1-phi)*tauc*C + taua*A - mu*R
    dD <- phi*tauc*C
    list(c(dS, dA, dC, dR, dD))}) 
}

parameters['rho']<-0.75
out3r_lowbeta <- deSolve::ode(y = istate_a, times = tps, func = SACRD_rho, parms = parameters) 

lines(out3r_lowbeta[,2],out3r_lowbeta[,3]+out3r_lowbeta[,4], col="violetred",lwd=3)




