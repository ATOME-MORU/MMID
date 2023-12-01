library(deSolve)
setwd('C:/mgh - mmid')

SEACRD<-function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    P <- (S+E+A+C+R)
    lam<-parameters['beta']*(rho*A+C+delta*E)/P
    
    dS <- mu*P-mu*S-lam*S   
    dE <- lam*S - v*E - mu*E
    dA <- (1-gamma)*v*E - taua*A - mu*A
    dC <- gamma*v*E - tauc*C - mu*C 
    dR <- (1-phi)*tauc*C + taua*A - mu*R
    dD <- phi*tauc*C
    
    list(c(dS, dE, dA, dC, dR, dD))}) 
}


initP<-999 
initE<-0
initA<-0
initC<-1
initR<-0 
initD<-0
initS<-initP-initE-initA-initC-initR
istate <- c(S = initS, E = initE, A = initA, C = initC, R = initR, D = initD)

parameters <- c(
  mu = (1/(50*52*7)),    
  beta = 0.25, 
  v = 1/4,
  gamma = 0.2,
  taua = 1/10,
  tauc = 1/10,
  phi=0.05,
  rho=1,
  delta=0.0
)

time_start <- 0
time_stop <- 600
deltat<- 0.1
tps <- seq(time_start , time_stop , by = deltat)

out <- deSolve::ode(y = istate, times = tps, func = SEACRD, parms = parameters) 
Result<-as.data.frame(out)

# plot(out)
# plot(Result)
plot(Result$S,Result$C+Result$A, ylim = c(0,300),xlim = c(0,1000))


parameters['v']<-1e10
out2 <- deSolve::ode(y = istate, times = tps, func = SEACRD, parms = parameters) 
Result2<-as.data.frame(out2)
lines(Result2$S,Result2$C+Result2$A, col = "lightcoral", lwd=3)


parameters['v']<-1/4
parameters['delta']<-0.5
out3 <- deSolve::ode(y = istate, times = tps, func = SEACRD, parms = parameters) 
Result3<-as.data.frame(out3)
lines(Result3$S,Result3$C+Result3$A, col = "plum", lwd=3)


plot(Result$R,Result$C+Result$A, ylim = c(0,300),xlim = c(0,1000))
lines(Result2$R,Result2$C+Result2$A, col = "lightcoral", lwd=3)
lines(Result3$R,Result3$C+Result3$A, col = "plum", lwd=3)


plot(Result$C+Result$A, ylim = c(0,300))
lines(Result2$C+Result2$A, col = "lightcoral", lwd=3)
lines(Result3$C+Result3$A, col = "plum", lwd=3)
sum(Result$C+Result$A)
sum(Result2$C+Result2$A)
sum(Result3$C+Result3$A)

v_vector<-seq(0.1,1,by=0.1)
delta_vector<-(seq(0.2,1,by=0.1))
result<-matrix(0,nrow = length(v_vector),length(delta_vector))
for (i in 1:length(v_vector)){
  for (j in 1:length(delta_vector)){
    parameters["v"]<-v_vector[i]
    parameters["delta"]<-delta_vector[j]
    
    out <- deSolve::ode(y = istate, times = tps, func = SEACRD, parms = parameters) 
    
    result[i,j]<-tail(out[,6],1)
  }
}
contour(result,xaxt = "n", yaxt = "n",xlab="v",ylab="delta",nlevels = 10, labcex=1.6,col = hcl.colors(12, "Temps"),lwd =2)
axis(1, at=1:length(v_vector)/length(v_vector), labels=v_vector)
axis(2, at=1:length(delta_vector)/length(delta_vector), labels=delta_vector)





