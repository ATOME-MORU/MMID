initS <- 100
initI <- 10
beta <- 0.5
gamma<-1/5
tau<-1/10
alpha<-1
no.of.timesteps <- 100
pop <- c(rep(0,initS),rep(1,initI))
sim.table <- as.data.frame(matrix(NA, no.of.timesteps, 5))
colnames(sim.table) <- c('Day','S','E','I','R')

for(t in 1:no.of.timesteps)
{
  #Simulate what happens at each timestep
  sim.table$Day[t] <- t
  sim.table$S[t] <- sum(pop==0)
  sim.table$E[t] <- sum(pop==1)
  sim.table$I[t] <- sum(pop==2)
  sim.table$R[t] <- sum(pop==3)
  lambda <- beta*sum(pop==2)/length(pop)
  
  random.lambda <- runif(length(pop))
  random.gamma <- runif(length(pop))
  random.tau <- runif(length(pop))
  random.alpha <- runif(length(pop))
  
  for(i in 1:length(pop)){
    #Simulate what happens to each individual
    
    ## INFECTION ##
    if(pop[i]==0){
      if(random.lambda[i]<lambda){
        pop[i] <- 1
      }
    } 
    
    ## LATENCY ##
    if(pop[i]==1){
      if(random.gamma[i]<gamma){ 
        pop[i] <- 2
      }
    }
    ## RECOVERY ##
    if(pop[i]==2){
      if(random.tau[i]<tau){
        pop[i] <- 3
      }
    }
    
    ## LOSS OF IMMUNITY ##
    if(pop[i]==3){
      if(random.alpha[i]<alpha){
        pop[i] <- 0
      }
    }
    
  }
}
sim.table

plot(sim.table$Day,sim.table$S,type='l',ylim=c(0,length(pop)),col='black',xlab="Weeks",ylab=("# People"))
lines(sim.table$Day,sim.table$I,col='red')
lines(sim.table$Day,sim.table$R,col='blue')
legend("topright", legend = c("S","I","R"), lwd = 1, col=c("black","red","blue"))
