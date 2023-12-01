
# Mathematical Modelling of Infectious Diseases
###################################################
## MODELLING INTERVENTIONS IN R PRACTICAL SESSION## 
###################################################

# Extend SEIR Model to include interventions
# Reduced risk
# vaccaination
# treatment 

## change to your WD 
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
                beta=1.75,               # basic reproduction number
                omega=(1/(10*52)),  # rate of loss of immunity = 1/(average duration of immunity)
                gamma=1/2,          # rate of movement from latent to infectious stage = 1/(average latent period)
                nui=1/7,            # rate of recovery = 1/(average duration of infection)
                report=1/7,         # proportion of all infections that are reported
                amp=0,              # relative amplitude of seasonal forcing
                phi=0,              # week of peak in seasonal forcing
                improv_treat = 1/2, # rate of treatment with drug
                gammatreat = 1/2,   # rate of recovery if using the drug
                
                week_interv = 5,   # weeks after first case when the intervention starts
                lat_cov_i = 0.24,   # intervention coverage of latrines
                lat_eff = 0.7,      # efficacy of latrines in reducing transmission
                vac_cov_i = 0.0,     # coverage of vaccine campaign
                prop_treat = 0.25,    # proportion treated
                infect_treat = 0.25, # reduced infectiousness while treated
                campaignweeks = 4   # time taken to reach target coverage
)

# MODEL INITIAL CONDITIONS
initP<-3700 # population size

initE<-1 # Exposed
initI<-0 # Infectious
initR<-0 # Immune
initT<-0 # Treated 
initS<-initP-initE-initI-initR+initT  # Susceptible (non-immune)

state <- c(S = initS, E=initE, I = initI,R = initR, TRT = initT)

# set up a function to solve the equations
HepE<-function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         
         # define variables
         P <- (S+E+I+R+TRT)
         seas<-1+amp*cos(2*pi*(t-phi)/52)
         R0<-beta*gamma/(muo+nui*(1-prop_treat)+prop_treat*improv_treat)*(gamma+muo)
         
         ## IMPLEMENT RISK REDUCING INTERVENTION AS A PULSE (LATRINES, BED NETS, E.G)
         lat_cov <- (t>=(week_start+week_interv))*lat_cov_i
         latrine <- (1-lat_eff*lat_cov)
         lam <- latrine*beta*seas*(I/P+infect_treat*TRT/P)
         
         # IMPLEMENT VECCINATION AS A RATE
         vac_rate <- (-log(1-vac_cov_i)/campaignweeks)
         vaccinate <- (t>=(week_start+week_interv))*(t<=week_start+week_interv+campaignweeks)*vac_rate
         
         # rate of change
         dS <- mui*P-muo*S-lam*S+omega*R-vaccinate*S
         dE <- -muo*E+lam*S-gamma*E
         dI <- -muo*I-(1-prop_treat)*nui*I-prop_treat*improv_treat*I+gamma*E
         dR <- -muo*R-omega*R+vaccinate*S+gammatreat*TRT+(1-prop_treat)*nui*I
         dTRT <- -muo*TRT-gammatreat*TRT+prop_treat*improv_treat*I
         
         # return the rate of change
         list(c(dS, dE, dI, dR, dTRT))
       }
  ) 
  
}


# run the model
out <- ode(y = state, times = times, func = HepE, parms = parameters)

# some more model outputs
# total population
pop<-out[,"S"]+out[,"E"]+out[,"I"]+out[,"R"]+out[,"TRT"]
# weekly incidence
inc <- parameters["report"]*parameters["gamma"]*out[,"E"]

time<-out[,"time"]

# COMPARING THE MODEL OUTPUT WITH DATA
# label axes 
plot(time,inc,type='l',lwd=3,main = "Predicted Hepatitis E Outbreak",xlab = "Time in weeks",ylab="New reported cases per week",ylim=c(0,max(data[,"Cases"],inc)))
# plot the data with the model output
points(data[,"week"],data[,"Cases"],pch=19,col='red')


### PLOT model variables
plot(time,out[,"S"],type='l',col="blue",lwd=3,main = "Model Variables",xlab = "Time in weeks",ylab="Number of people")
lines(time,out[,"E"],type='l',col="yellow2",lwd=3)
lines(time,out[,"I"],type='l',col="red2",lwd=3)
lines(time,out[,"R"],type='l',col="green3",lwd=3)
lines(time,out[,"TRT"],type='l',col="orange1",lwd=3)
legend(65, 2000, legend=c("S","E","I","R","TRT"),
       col=c("blue", "yellow2","red2","green3","orange1"), lty=1,lwd=2, cex=1.3)

## PLOT TOTAL INCIDENCE
sum_inc<-cumsum(inc)
plot(time,sum_inc,type='l',lwd=3)

