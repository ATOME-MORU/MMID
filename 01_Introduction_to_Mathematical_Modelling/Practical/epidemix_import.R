library(deSolve)
setwd('C:/mgh - mmid')

app<-read.csv('SIR.csv', header = T)
app

plot(app$S, app$Is)
plot(app)
plot(app[,c(2,5,6)])


SIR_model<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+I+R)
    lam<-parameters['beta']*I/P
    
    dS <- mu*P-mu*S-lam*S         
    dI <- lam*S-tau*I-mu*I 
    dR <- tau*I-mu*R
    list(c(dS, dI, dR))}) 
}


initP<-9999
initI<-1 
initR<-0 
initS<-initP-initI-initR
istate <- c(S = initS, I = initI, R = initR)

parameters <- c(
  mu=(1/(50*52*7)),    
  beta=0.4,              
  tau=1/10            
)

time_start <- 0
time_stop <- 365*2 
deltat<- 0.1
tps <- seq(time_start , time_stop , by = deltat)

out <- deSolve::ode(y = istate, times = tps, func = SIR_model, parms = parameters) 
plot(out[,1],out[,3])
lines(app$X,app$Is, col="red")

out2 <- deSolve::ode(y = istate, times = app$X, func = SIR_model, parms = parameters) 
plot(out2[,1],out2[,3])
lines(app$X,app$Is, col="red")



SIRS_model<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+I+R)
    lam<-parameters['beta']*I/P
    
    dS <- mu*P-mu*S-lam*S+alpha*R         
    dI <- lam*S-tau*I-mu*I 
    dR <- tau*I-mu*R-alpha*R
    list(c(dS, dI, dR))}) 
}

parameters['alpha']<-1/150 
out <- deSolve::ode(y = istate, times = tps, func = SIRS_model, parms = parameters) 
plot(out[,1],out[,3])

parameters['beta']<-parameters['beta']*5
out_high_beta <- deSolve::ode(y = istate, times = tps, func = SIRS_model, parms = parameters) 
lines(out_high_beta[,1],out_high_beta[,3],col="red")

parameters['beta']<-parameters['beta']/10
out_low_beta <- deSolve::ode(y = istate, times = tps, func = SIRS_model, parms = parameters) 
lines(out_low_beta[,1],out_low_beta[,3],col="blue")

parameters['beta']<-0.4
parameters['alpha']<-parameters['alpha']*5
out_high_alpha <- deSolve::ode(y = istate, times = tps, func = SIRS_model, parms = parameters) 
lines(out_high_alpha[,1],out_high_alpha[,3],col="magenta")

parameters['alpha']<-parameters['alpha']/10
out_low_alpha <- deSolve::ode(y = istate, times = tps, func = SIRS_model, parms = parameters) 
lines(out_low_alpha[,1],out_low_alpha[,3],col="green3")


parameters['alpha']<-1/150 
out_sir <- deSolve::ode(y = istate, times = tps, func = SIR_model, parms = parameters) 
out_sirs <- deSolve::ode(y = istate, times = tps, func = SIRS_model, parms = parameters) 

plot(out_sir[,2],out_sir[,3],xlab = "S", ylab = "I")
lines(out_sirs[,2],out_sirs[,3],col="red")
lines(rep(2500,4000),1:4000)

parameters['alpha']<-1/50 
out_sirs2 <- deSolve::ode(y = istate, times = tps, func = SIRS_model, parms = parameters) 
lines(out_sirs2[,2],out_sirs2[,3],col="green3")


