
###################################################
## MODELLING INTERVENTIONS IN R PRACTICAL SESSION## 
###################################################

setwd("C:/mgh - mmid")

library(deSolve)
data <- read.csv('IHTM_hepEdata.csv')

# Set the start and end time for the model simulation
week_start <- 38
week_stop <- 73
times <- seq(week_start, week_stop, by = (1/7))

#MODEL PARAMETERS
parameters <- c(mui=(1/(50*52)),    # birth
                muo=(1/(50*52)),    # death
                R0=6.6,               # basic reproduction number
                omega=(1/(10*52)),  # rate of loss of immunity = 1/(average duration of immunity)
                gamma=1/2,          # rate of movement from latent to infectious stage = 1/(average latent period)
                nui=1/4,            # rate of recovery = 1/(average duration of infection)
                report=1/7,         # proportion of all infections that are reported
                amp=0,              # relative amplitude of seasonal forcing
                phi=0,              # week of peak in seasonal forcing
                
                week_interv = 15,   # weeks after first case when the intervention starts
                lat_cov_i = 0.24,   # intervention coverage of latrines
                lat_eff = 0.7       # efficacy of latrines in reducing transmission
)

# MODEL INITIAL CONDITIONS
initP<-3700 # population size

initE<-1 # Exposed
initI<-0 # Infectious
initR<-0 # Immune
initS<-initP-initE-initI-initR # Susceptible (non-immune)

state <- c(S = initS, E=initE, I = initI,R = initR)

# set up a function to solve the equations
HepE<-function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         
         # define variables
         P <- (S+E+I+R)
         seas<-1+amp*cos(2*pi*(t-phi)/52)
         beta<-R0*(muo+nui)*(gamma+muo)/gamma
         
         lat_cov <- (t>=(week_start+week_interv))*lat_cov_i
         latrine <- (1-lat_eff*lat_cov)
        
         lam <- latrine*beta*seas*I/P
         
         # rate of change
         dS <-  mui*P-muo*S-lam*S+omega*R
         dE <- -muo*E+lam*S-gamma*E
         dI <- -muo*I+gamma*E-nui*I
         dR <- -muo*R+nui*I-omega*R
         
         # return the rate of change
         list(c(dS, dE, dI, dR))
       }
  ) 
  
}


# run the model
out <- ode(y = state, times = times, func = HepE, parms = parameters)

# some more model outputs
# total population
pop<-out[,"S"]+out[,"E"]+out[,"I"]+out[,"R"]
# weekly incidence
inc <- parameters["report"]*parameters["gamma"]*out[,"E"]

time<-out[,"time"]

# COMPARING THE MODEL OUTPUT WITH DATA
# label axes 
plot(time,inc,type='l',lwd=3,main = "Predicted Hepatitis E Outbreak",xlab = "Time in weeks",ylab="New reported cases per week",ylim=c(0,max(data[,"Cases"],inc)))
# plot the data with the model output
points(data[,"week"],data[,"Cases"],pch=19,col='red')


