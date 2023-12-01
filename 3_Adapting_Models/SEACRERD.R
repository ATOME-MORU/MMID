library(deSolve)
setwd('C:/mgh - mmid')

SEACRERD<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+E+A+C+R+ER)
    lam<-parameters['beta']*(rho*A+C+delta*E+delta_r*ER)/P
    
    dS <- mu*P-mu*S-lam*S + alpha*R  
    dE <- lam*S - v*E - mu*E
    dA <- (1-gamma)*v*E - taua*A - mu*A + v*ER
    dC <- gamma*v*E - tauc*C - mu*C 
    dR <- (1-phi)*tauc*C + taua*A - mu*R - sigma*lam*R - alpha*R
    dER <- sigma*lam*R - v*ER - mu*ER
    dD <- phi*tauc*C
    
    list(c(dS, dE, dA, dC, dR, dER, dD))}) 
}


initP<-999 
initE<-0
initA<-0
initC<-1
initR<-0 
initER<-0
initD<-0
initS<-initP-initE-initA-initC-initR-initER
istate <- c(S = initS, E = initE, A = initA, C = initC, R = initR, ER = initER, D = initD)

parameters <- c(
  mu = (1/(50*52*7)),    
  beta = 0.25, 
  v = 1/4,
  gamma = 0.2,
  taua = 1/10,
  tauc = 1/10,
  phi = 0.05,
  rho = 1,
  sigma = 0.5,
  alpha = 1/180,
  delta_r = 0.0,
  delta=0.0
)

time_start <- 0
time_stop <- 600
deltat<- 0.1
tps <- seq(time_start , time_stop , by = deltat)

out <- deSolve::ode(y = istate, times = tps, func = SEACRERD, parms = parameters) 
Result<-as.data.frame(out)

plot(out)
# plot(Result)
plot(Result$S,Result$C+Result$A, ylim = c(0,300),xlim = c(0,1000))


beta_vector<-seq(0.1,3,by=0.2)
alpha_vector<-(seq(0.003,0.01,by=0.001))
result<-matrix(0,nrow = length(beta_vector),length(alpha_vector))
for (i in 1:length(beta_vector)){
  for (j in 1:length(alpha_vector)){
    parameters["beta"]<-beta_vector[i]
    parameters["alpha"]<-alpha_vector[j]
    
    out <- deSolve::ode(y = istate, times = tps, func = SEACRERD, parms = parameters) 
    
    result[i,j]<-tail(out[,5],1)
  }
}
contour(result,xaxt = "n", yaxt = "n",xlab="beta",ylab="alpha",nlevels = 10, labcex=1.6,col = hcl.colors(12, "Temps"),lwd =2)
axis(1, at=1:length(beta_vector)/length(beta_vector), labels=beta_vector)
axis(2, at=1:length(alpha_vector)/length(alpha_vector), labels=alpha_vector)




