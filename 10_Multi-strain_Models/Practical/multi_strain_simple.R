library(deSolve)
setwd('C:/mgh - mmid')

SIR_2strains<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- S+I1+I2+R1+R2+I12+I21+R12
    lam1<-beta*I1/P
    lam2<-beta*I2/P
    
    dS <- mu*P - mu*S- lam1*S - lam2*S + alpha1*R1 + alpha2*R2 
    dI1 <- lam1*S - tau1*I1 - mu*I1 
    dR1 <- tau1*I1 - mu*R1 - sigma*lam2*R1 - alpha1*R1 + alpha2*R12
    dI2 <- lam2*S - tau2*I2 - mu*I2 
    dR2 <- tau2*I2 - mu*R2 - sigma*lam1*R2 - alpha2*R2 + alpha1*R12
    dI12<- sigma*lam2*R1 - tau2*I12 - mu*I12
    dI21<- sigma*lam1*R2 - tau1*I21 - mu*I21
    dR12<- tau2*I12+tau1*I21 - mu*R12 - alpha1*R12 - alpha2*R12
    
    list(c(dS, dI1, dR1, dI2, dR2, dI12, dI21, R12))}) 
}


initP<-999 
initI1<-1
initR1<-0
initR2<-0
initI2<-1
initI12<-0
initI21<-0
initR12<-0
initS<-initP-initI1-initR1-initR2-initI2-initI12-initI21-initR12
istate <- c(S = initS, I1 = initI1, R1 = initR1, 
            I2 = initI2, R2 = initR2, I12 = initI12, I21=initI21,
            R12 = initR12)

parameters <- c(
  mu = (1/(50*52*7)),    
  beta = 0.5,
  tau1 = 1/10,
  tau2 = 1/10,
  sigma= 0.2,
  alpha1= 1/300,
  alpha2= 1/365
)

time_start <- 0
time_stop <- 1000
deltat<- 0.01
tps <- seq(time_start , time_stop , by = deltat)


out <- deSolve::ode(y = istate, times = tps, func = SIR_2strains, parms = parameters) 
plot(out[,'I1']+out[,'I21'])
lines(out[,'I2']+out[,'I12'], col = "cyan")


istate['I1']<-10
out2 <- deSolve::ode(y = istate, times = tps, func = SIR_2strains, parms = parameters) 
plot(out2[,'I1']+out2[,'I21'],ylim=c(0,350))
lines(out2[,'I2']+out2[,'I12'], col = "cyan")


parameters['tau2']<-1/20
out3 <- deSolve::ode(y = istate, times = tps, func = SIR_2strains, parms = parameters) 
lines(out3[,'I2']+out3[,'I12'], col = "magenta")


parameters['sigma']<-0
out4 <- deSolve::ode(y = istate, times = tps, func = SIR_2strains, parms = parameters) 
plot(out4[,'I1']+out4[,'I21'])
lines(out4[,'I2']+out4[,'I12'], col = "purple")

plot(out2[,'I1']+out2[,'I21'])
lines(out4[,'I1']+out4[,'I21'],col="red")
