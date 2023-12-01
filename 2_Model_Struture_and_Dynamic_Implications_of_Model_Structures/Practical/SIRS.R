library(deSolve)
setwd('C:/mgh - mmid')

SIRS_model<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+I+R)
    lam<-parameters['beta']*I/P
    
    dS <- mu*P-mu*S-lam*S+alpha*R         
    dI <- lam*S-tau*I-mu*I 
    dR <- tau*I-mu*R-alpha*R
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
  tau=1/10,
  alpha=1/150
)

R0<-parameters['beta']/
  (parameters['tau']+parameters['mu'])

time_start <- 0
time_stop <- 365 
deltat<- 0.1
tps <- seq(time_start , time_stop , by = deltat)

out <- deSolve::ode(y = istate, times = tps, func = SIRS_model, parms = parameters) 
plot(out)
tail(out[,3],1)


beta_vector<-seq(0.5,1.5,by=0.25)*parameters["beta"]
alpha_vector<-round(seq(0.5,1.5,by=0.25)*parameters["alpha"],4)
result<-matrix(0,nrow = length(beta_vector),length(alpha_vector))
for (i in 1:length(beta_vector)){
  for (j in 1:length(alpha_vector)){
    parameters["beta"]<-beta_vector[i]
    parameters["alpha"]<-alpha_vector[j]
    
    out <- deSolve::ode(y = istate, times = tps, func = SIRS_model, parms = parameters) 
    
    result[i,j]<-tail(out[,3],1)
  }
}

contour(result,xaxt = "n", yaxt = "n",ylab="alpha",xlab="beta",nlevels = 10,col = hcl.colors(14, "Temps"),lwd =2)
axis(1, at=1:length(beta_vector)/length(beta_vector), labels=beta_vector)
axis(2, at=1:length(alpha_vector)/length(alpha_vector), labels=alpha_vector)

library(plotly)
fig <- plot_ly(x=beta_vector,y=alpha_vector,z = ~result, type = "contour") %>%
  colorbar(title='Infectious people at equilibrium', titleside="right",x=0.91, ticks = "inside") %>%
  layout(xaxis = list(title = 'beta'),
         yaxis = list(title = 'alpha'))
fig



#############################################################################################################################
######################################################################################################################################


#                                       R0  as an input  parameter

###############################################################################################################################
###############################################################################################################################

SIRS_model2<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+I+R)
    lam<-(R0*tau)*I/P
    
    dS <- mu*P-mu*S-lam*S+alpha*R         
    dI <- lam*S-tau*I-mu*I 
    dR <- tau*I-mu*R-alpha*R
    list(c(dS, dI, dR))}) 
}


initP<-999 
initI<-1 
initR<-0 
initS<-initP-initI-initR
istate <- c(S = initS, I = initI, R = initR)

parameters2 <- c(
  mu=(1/(50*52*7)),    
  R0=2,              
  tau=1/10,
  alpha=1/150
)

time_start <- 0
time_stop <- 600 
deltat<- 0.1
tps <- seq(time_start , time_stop , by = deltat)

out2 <- deSolve::ode(y = istate, times = tps, func = SIRS_model2, parms = parameters2) 
plot(out2)
tail(out2[,3],1)


r0_vector<-seq(0.1,2,by=0.2)*parameters2["R0"]
result2<-rep(0,nrow = length(r0_vector))
for (i in 1:length(r0_vector)){
    parameters2["R0"]<-r0_vector[i]
    out2 <- deSolve::ode(y = istate, times = tps, func = SIRS_model2, parms = parameters2) 
    result2[i]<-tail(out2[,3],1)
}
plot(r0_vector,result2,type = "l", ylab = "Infectious people at equilibrium", xlab = "Basic reproduction number, R0")

