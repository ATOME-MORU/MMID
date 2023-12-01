library(deSolve)

##############################################################################

#            Competing death and recovery rates

#############################################################################

SIR_model<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+I+R)
    lam<-parameters['beta']*I/P
    
    dS <- mu*P-mu*S-lam*S         
    dI <- lam*S-tau*I-delta*I-mu*I 
    dR <- tau*I-mu*R
    dD <- delta*I
    list(c(dS, dI, dR, dD))}) 
}


initP<-999 
initI<-1 
initR<-0 
initD<-0
initS<-initP-initI-initR
istate <- c(S = initS, I = initI, R = initR, D = initD)

parameters <- c(
  mu=(1/(50*52*7)),    
  beta=0.25,              
  tau=1/10,
  delta = 0.25*1/10
)

time_start <- 0
time_stop <- 600
deltat<- 0.1
tps <- seq(time_start , time_stop , by = deltat)

out <- deSolve::ode(y = istate, times = tps, func = SIR_model, parms = parameters) 
Result<-as.data.frame(out)

# plot(Result)
plot(Result$S,Result$I)
plot(out)


##############################################################################

#            Mutually exclusive death and recovery rates

###############################################################################

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

parameters2 <- c(
  mu = (1/(50*52*7)),    
  beta = 0.25,              
  tau = 1/10,
  delta = 1/10,
  phi=0.25
)

R0<-parameters2['beta']/
  (parameters2['tau']+parameters2['mu'])
R0

out2 <- deSolve::ode(y = istate, times = tps, func = SIR_model2, parms = parameters2) 
Result2<-as.data.frame(out2)

plot(out2)
# plot(Result)
plot(Result2$S,Result2$I)

plot(Result$time,Result$D, ylim = c(0,400))
lines(Result2$time,Result2$D,col="red")


#################   Which one gives you the desired disease fatality rate?  ###########

tail(Result$D)/(tail(Result$D)+tail(Result$R))
tail(Result2$D)/(tail(Result2$D)+tail(Result2$R))



#################         Imported infections         ########################################    

SIR_model2_import<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P <- (S+I+R)
    lam<-parameters['beta']*(I+imp)/P
    
    dS <- mu*P-mu*S-lam*S         
    dI <- lam*S-(1-phi)*tau*I-delta*phi*I-mu*I 
    dR <- (1-phi)*tau*I-mu*R
    dD <- delta*phi*I
    list(c(dS, dI, dR, dD))}) 
}

parameters2['imp']<-10

out2_import <- deSolve::ode(y = istate, times = tps, func = SIR_model2_import, parms = parameters2) 
Result2_import<-as.data.frame(out2_import)

plot(Result2$S,Result2$I, ylim = c(0, 300))
lines (Result2_import$S,Result2_import$I, col = "magenta")



parameters2['beta']<-0.09
import_vector<-seq(0,10,by=1)
result_import<-matrix(0,nrow = 1,ncol = 3)
colnames(result_import)<-c("time","I","imp")
result_import<-as.data.frame(result_import)
for (i in 1:length(import_vector)){
  parameters2["imp"]<-import_vector[i]
  out_aux <- deSolve::ode(y = istate, times = tps, func = SIR_model2_import, parms = parameters2) 
  aux_mat<-cbind(out_aux[,c(1,3)],rep(import_vector[i],length(tps)))
  colnames(aux_mat)<-c("time","I","imp")
  result_import<-rbind(result_import,aux_mat)
}

library(ggplot2)
g<-ggplot(result_import)
g<-g+geom_line(aes(x=time,y=I,colour=as.factor(imp)))
g

