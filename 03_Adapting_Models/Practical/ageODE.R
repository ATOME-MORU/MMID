###########################################################################
## AGE-DEPENDANT SEIRS MODEL WITH 5-YEAR AGE CLASSES  with UN demog data ##
###########################################################################

require("deSolve")
setwd("C:/mgh2022_23")

###########################################################################
# DEMOGRAPHIC DATA
# https://population.un.org/wpp/Download/Standard/Population/

# population structure in 2020
popstruc <- read.csv("pop_structure_X.csv", header = TRUE)
# convert from 1000s to total numbers
popstruc[, 2] <- popstruc[, 2] * 1000
A <- length(popstruc[, 2])

# births by age of mother
popbirth <- read.csv("pop_birth_X.csv", header = TRUE)
# convert from 1000s per 5 year period to per person per day
popbirth[, 2] <- 1000 * popbirth[, 2] / (5 * popstruc[, 2] * 365.25)

#natural mortality per person per year
popmort <- read.csv("pop_mort_X.csv", header = TRUE)
# convert from 1000s per 5 year period to per person per day
popmort[, 2] <- 1000 * popmort[, 2] / (5 * popstruc[, 2] * 365.25)
mort <- popmort[, 2]



###########################################################################
# CONTACT DATA
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697

c_home <- as.matrix(read.csv("contact_home_X.csv", header = FALSE))
c_school <- as.matrix(read.csv("contact_school_X.csv", header = FALSE))
c_work <- as.matrix(read.csv("contact_work_X.csv", header = FALSE))
c_other <- as.matrix(read.csv("contact_other_X.csv", header = FALSE))
nce <- A - length(c_home[1, ])

# filling in 4 higher age groups 75-80, 80-85, 85-90, 95-100, 100+
contact_home <- matrix(0, nrow = A, ncol = A)
contact_school <- matrix(0, nrow = A, ncol = A)
contact_work <- matrix(0, nrow = A, ncol = A)
contact_other <- matrix(0, nrow = A, ncol = A)

for (i in 1:(A - nce)){
  for (j in 1:(A - nce)){
    contact_home[i, j] <- c_home[i, j]
    contact_school[i, j] <- c_school[i, j]
    contact_work[i, j] <- c_work[i, j]
    contact_other[i, j] <- c_other[i, j]
  }
}

for (i in (A + 1 - nce):A){
  for (j in 1:(A - nce)){
    contact_home[i, j] <- c_home[(A - nce), j]
    contact_school[i, j] <- c_school[(A - nce), j]
    contact_work[i, j] <- c_work[(A - nce), j]
    contact_other[i, j] <- c_other[(A - nce), j]
  }
}
for (i in 1:(A - nce)){
  for (j in (A + 1 - nce):A){
    contact_home[i, j] <- c_home[i, (A - nce)]
    contact_school[i, j] <- c_school[i, (A - nce)]
    contact_work[i, j] <- c_work[i, (A - nce)]
    contact_other[i, j] <- c_other[i, (A - nce)]
  }
}
for (i in (A + 1 - nce):A){
  for (j in (A + 1 - nce):A){
    contact_home[i, j] <- c_home[(A - nce),(A - nce)]
    contact_school[i, j] <- c_school[(A - nce),(A - nce)]
    contact_work[i, j] <- c_work[(A - nce),(A - nce)]
    contact_other[i, j] <- c_other[(A - nce),(A - nce)]
  }
}

###########################################################################
# average contacts per day from POLYMOD matrices 
c <- sum((contact_home + contact_other + contact_school + contact_work) %*%
      (popstruc[, 2] / sum(popstruc[, 2])))


#####  IHR and CFR from data
# of confirmed cases requiring ICU
csr <- read.csv("covidagesevere_X.csv", header = TRUE)
# of confirmed cases who die
cfr <- read.csv("covidagemort_X.csv", header = TRUE)
# csr-cfr >= 0 for all age classes - if not model won't work
for (i in 1 : A) {
  cfr[i, 2] <- min(csr[i, 2], cfr[i, 2])
}


###########################################################################
# per year ageing matrix
dd <- seq(1:A) / seq(1:A)
ageing <- t(diff(diag(dd), lag = 1) / (5 * 365.25))
ageing <- cbind(ageing, 0 * seq(1:A)) # no ageing from last compartment
###########################################################################


initP <- sum(popstruc[, 2]) # population size
ageindcase <- 20 # age of index case (years)
aci <- floor((ageindcase / 5) + 2) # age class of index case

startdate <- as.Date("2020-01-31")
stopdate <- as.Date("2020-12-31")
day_start <- as.numeric(startdate-startdate)
day_stop <- as.numeric(stopdate-startdate)
times <- seq(day_start, day_stop)


#MODEL PARAMETERS
parameters <- c(
  p = 0.047,                # probabilty of infection given a contact 
  rho = 0.2,                # relative infectiousness of incubation phase
  omega = (1 / (100 * 365)),# rate of loss of immunity = 1/(average duration of immunity) # nolint
  gamma = 1 / 4.5,          # rate of movement from incubation to infectious stage = 1/(average incubation period)
  nui = 1 / 5,              # rate of recovery = 1/(average duration of symptomatic infection phase)
  report = 1 / 8,           # proportion of all infections that are reported
  ratem = 1 / 7,            # 1/ time to death for fatal
  nus = 1 / 7,              # rate if recovery for non-fatal severe
  rhos = 0.1               # relative infectiousness*contacts of severe ill patients 
)


###########################################################################
# Define the indices for each variable
Sindex <- 1 : A
Eindex <- (A + 1) : (2 * A)
Iindex <- (2 * A + 1) : (3 * A)
Rindex <- (3 * A + 1) : (4 * A)
Hindex <- (4 * A + 1) : (5 * A)
Mindex <- (5 * A + 1) : (6 * A)
Cindex <- (6 * A + 1) : (7 * A)
CMindex <- (7 * A + 1) : (8 * A)

###########################################################################
# MODEL INITIAL CONDITIONS
initI <- 0 * popstruc[, 2]# Infected and symptomatic
initE <- 0 * popstruc[, 2] # Incubating
initE[aci] <- 1 # place random index case in E compartment
initR <- 0 * popstruc[, 2] # Immune
initH <- 0 * popstruc[, 2] # hospitalised 
initM <- 0 * popstruc[, 2] # died 
initC <- 0 * popstruc[, 2] # Cumulative cases (true)
initCM <- 0 * popstruc[, 2] # Cumulative deaths (true)
initS <- popstruc[, 2] - initE - initI -
          initR - initH - initM # Susceptible (non-immune)

# initial conditions for the main solution vector
Y <- c(initS, initI, initE, initR,
       initH, initM, initC, initCM)

# set up a function to solve the equations
covid <- function(t, Y, parameters) {
  with(as.list(c(Y, parameters)),
    {
      S <- Y[Sindex]
      E <- Y[Eindex]
      I <- Y[Iindex]
      R <- Y[Rindex]
      H <- Y[Hindex]
      M <- Y[Mindex]
      C <- Y[Cindex]
      CM <- Y[CMindex]
      # print(t)
      contacts <- contact_home + contact_other + contact_school + contact_work
      # print(dim(contacts))
      P <- (S + E + I + R + H + M)
      b1 <- sum(popbirth[, 2] * P)
      birth <- 0 * popbirth[, 2]
      birth[1] <- b1
         
      lam <- p * contacts %*% ((rho * E + I + rhos * H) / P)
      # print(dim(lam))
      dSdt <- -S * lam + omega * R + ageing %*% S - mort * S + birth
      dEdt <- S * lam - gamma * E + ageing %*% E - mort * E
      dIdt <- gamma * (1 - csr[, 2]) * E - nui * I + ageing %*% I - mort * I
      dRdt <- nui * I - omega * R + ageing %*% R - mort * R + nus * H

      dHdt <- ageing %*% H + gamma * (csr[, 2] - cfr[, 2]) * E -
               mort * H - nus * H
      dMdt <- ageing %*% M + gamma * cfr[, 2] * E -  ratem * M - mort * M
         
      dCdt <- report * gamma * E
      dCMdt <- ratem * M + mort * M
         
      # return the rate of change
      list(c(dSdt, dEdt, dIdt, dRdt, dHdt, dMdt, dCdt, dCMdt))
    }
  )
}

# run the model
out <- ode(y = Y, times = times, func = covid, parms = parameters)

# total population
pop <- out[, (Sindex + 1)] + out[, (Eindex + 1)] + out[, (Iindex + 1)] +
       out[, (Rindex + 1)] + out[, (Hindex + 1)]
tpop <- rowSums(pop)
time <- as.Date(out[, 1] + startdate)
# daily incidence
inc <- parameters["report"] * parameters["gamma"] * out[, (Eindex + 1)]
dailyinc <- rowSums(inc)
hospreq <- rowSums(out[, (Hindex + 1)]) # requirement for hospital beds
cmortality <- rowSums(out[, (CMindex + 1)]) # cumulative mortality

## plot incidence over time
plot(rowSums(inc))

# plot hospital occupancy against the number of hospital beds available
plot(hospreq)

#### age groups
deaths <- (out[, CMindex + 1])
colnames(deaths) <- c("0-4", "5-9", "10-14", "15-19", "20-24",
            "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
            "55-59", "60-64", "65-69", "70-74", "75-79", "80-84",
            "85-89", "90-94", "95-99", "100+")
deaths_prop <- round(tail(deaths, 1), 0) / sum(round(tail(deaths, 1), 0))
barplot(deaths_prop)


# output doubling time over time first 7 days
dd <- 7
doub0 <- log(2) * dd / (log(inc[2 + dd] / inc[2]))
doub0


# run the model for the intervention

# implement social distancing measures at specific times
parameters_intervention <- parameters
parameters_intervention["social_distancing_start"] <- 20
parameters_intervention["social_distancing_end"] <- 80
parameters_intervention["social_distancing_effect"] <- 0.15
parameters_intervention["shielding_start"] <- 1e100
parameters_intervention["shielding_end"] <- 80
shielding_effect <- seq(0, 100, by = 5)

# set up a function to solve the equations
covid_interventions <- function(t, Y, parameters) {
  with(as.list(c(Y, parameters)),
    {
      S <- Y[Sindex]
      E <- Y[Eindex]
      I <- Y[Iindex]
      R <- Y[Rindex]
      H <- Y[Hindex]
      M <- Y[Mindex]
      C <- Y[Cindex]
      CM <- Y[CMindex]
      # print(t)
      contacts <- contact_home + contact_other + contact_school + contact_work

      P <- (S + E + I + R + H + M)
      b1 <- sum(popbirth[, 2] * P)
      birth <- 0 * popbirth[, 2]
      birth[1] <- b1

      social_distancing <- ifelse(t > social_distancing_start &
        t < social_distancing_end,
        social_distancing_effect, 1)

      
      shielding <- if (t > shielding_start &
        t < shielding_end) {
          (1 - shielding_effect)
        }else {
           rep(1, A)
        }

      # print(shielding_effect)
      # print(paste("shield: ",shielding))
      # print(dim(shielding))
      # print(paste("t: ", t, "sd: ", social_distancing))

      lam <- p * social_distancing *
        contacts %*% (shielding * (rho * E + I + rhos * H) / P)
         
      dSdt <- -S * lam + omega * R + ageing %*% S - mort * S + birth
      dEdt <- S * lam - gamma * E + ageing %*% E - mort * E
      dIdt <- gamma * (1 - csr[, 2]) * E - nui * I + ageing %*% I - mort * I
      dRdt <- nui * I - omega * R + ageing %*% R - mort * R + nus * H

      dHdt <- ageing %*% H - mort * H + gamma * (csr[, 2] - cfr[, 2]) * E -
               mort * H - nus * H
      dMdt <- ageing %*% M + gamma * cfr[, 2] * E -  ratem * M - mort * M
         
      dCdt <- report * gamma * E
      dCMdt <- ratem * M + mort * M
         
      # return the rate of change
      list(c(dSdt, dEdt, dIdt, dRdt, dHdt, dMdt, dCdt, dCMdt))
    }
  )
}

# run the model
out_interv <- ode(y = Y, times = times, func = covid_interventions,
  parms = parameters_intervention)

# total population
pop_interv <- out_interv[, (Sindex + 1)] + out_interv[, (Eindex + 1)] + out_interv[, (Iindex + 1)] +
       out_interv[, (Rindex + 1)] + out_interv[, (Hindex + 1)]
tpop_interv <- rowSums(pop_interv)
time <- as.Date(out_interv[, 1] + startdate)
# daily incidence
inc_interv <- parameters_intervention["report"] * parameters_intervention["gamma"] *
  out_interv[, (Eindex + 1)]
dailyinc_interv <- rowSums(inc_interv)
hospreq_interv <- rowSums(out_interv[, (Hindex + 1)])
cmortality_interv <- rowSums(out_interv[, (CMindex + 1)]) # cumulative mortality

## plot incidence over time
plot(rowSums(inc_interv), type = "l", col = "blue")
lines(rowSums(inc), type = "l", col = "red")

# plot hospital occupancy against the number of hospital beds available
plot(hospreq_interv)

#### age groups
deaths_interv <- (out_interv[, CMindex + 1])
colnames(deaths_interv) <- c("0-4", "5-9", "10-14", "15-19", "20-24",
            "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
            "55-59", "60-64", "65-69", "70-74", "75-79", "80-84",
            "85-89", "90-94", "95-99", "100+")
deaths_prop_interv <- round(tail(deaths_interv, 1), 0) /
  sum(round(tail(deaths_interv, 1), 0))
barplot(deaths_prop_interv)

death_data <- data.frame("Age" = 1:21, "Control"=as.vector(round(tail(deaths, 1), 0)), "Intervention" = as.vector(round(tail(deaths_interv, 1), 0)))
library(reshape2)
df_long <- melt(death_data, id.var = "Age")
library(ggplot2)
ggplot(df_long, aes(x = Age, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") 



## compare with shielding
parameters_intervention["shielding_start"] <- 5
parameters_intervention["shielding_end"] <- 800
shielding_effect <- seq(0, 1, by = .05)

out_interv_shield <- ode(y = Y, times = times, func = covid_interventions,
  parms = parameters_intervention)
deaths_interv_shielding <- (out_interv_shield[, CMindex + 1])


death_data <- data.frame("Age" = 1:21, "Control" = as.vector(round(tail(deaths, 1), 0)), 
  "Lockdown" = as.vector(round(tail(deaths_interv, 1), 0)),
  "Shielding" = as.vector(round(tail(deaths_interv_shielding , 1), 0)))
library(reshape2)
df_long <- melt(death_data, id.var = "Age")
library(ggplot2)
ggplot(df_long, aes(x = Age, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge")


plot(rowSums(deaths), type = "l", col = "red")
lines(rowSums(deaths_interv), type = "l", col = "blue")
lines(rowSums(deaths_interv_shielding), type = "l", col = "seagreen")
legend("topright", legend = c("Control", "Lockdown", "Shielding"), lwd = 1,
    col = c("red", "blue", "seagreen"))
