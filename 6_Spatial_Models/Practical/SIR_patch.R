library(deSolve)
setwd("C:/mgh/spatial")

SIR_patch<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    P1 <- S1 + I1 + R1
    P2 <- S2 + I2 + R2

    if (t >= tbit){beta1 <- beta1 * tbeff}

    lam1 <- beta1 * I1 / P1
    lam2 <- beta2 * I2 / P2
    
    dS1 <- mu*P1 - mu*S1 - lam1*S1 - f12*S1 + f21*S2 
    dI1 <- lam1*S1 - tau*I1 - mu*I1 - f12*I1 + f21*I2
    dR1 <- tau*I1 - mu*R1 - f12*R1 + f21*R2
    dS2 <- mu*P2 - mu*S2 - lam2*S2 + f12*S1 - f21*S2  
    dI2 <- lam2*S2 - tau*I2 - mu*I2 + f12*I1 - f21*I2
    dR2 <- tau*I2 - mu*R2 + f12*R1 - f21*R2

    N <- P1 + P2
    count <- 1
    list(c(dS1, dI1, dR1, dS2, dI2, dR2, N, count))})
}


initP1 <- 1000
initI1 <- 1
initR1 <- 0
initS1 <- initP1 - initI1 - initR1

initP2 <- 1000
initR2 <- 0
initI2 <- 0
initS2 <- initP2 - initI2 - initR2

istate <- c(S1 = initS1, I1 = initI1, R1 = initR1, 
            S2 = initS2, I2 = initI2, R2 = initR2, 
            N = initS1 + initS1 + initR1 + initS2 + initS2 + initR2,
            cout = 0)

parameters <- c(
  mu = (1 / (50 * 52 * 7)),
  beta1 = 0.75,
  beta2 = 0.45,
  tau = 1 / 10,
  f12 = 0.001,
  f21 = 0.001,
  tbit = 1e100
)

time_start <- 0
time_stop <- 100
deltat<- 0.001
tps <- seq(time_start , time_stop , by = deltat)
tt <- seq(1, length(tps), by = 1/deltat)

out <- deSolve::ode(y = istate, times = tps, func = SIR_patch, parms = parameters) 
plot(out)
ww <- which(out[, 1] %% 1 == 0)
sum(out[ww, 2:4]) + sum(out[ww, 5:7])
tail(out[, 8], 1)
tail(out[, 9], 1)


parameters["f12"] <- 0.0
out_no_migration <- deSolve::ode(y = istate, times = tps, func = SIR_patch, parms = parameters) 
plot(out_no_migration)

parameters["f12"] <- 0.001
parameters["tbit"] <- 10
parameters["tbeff"] <- 0.2
out_intervention <- deSolve::ode(y = istate, times = tps, func = SIR_patch, parms = parameters) 
par(mfrow=c(2,1))
plot(out[,3], type = "l", col = "orange3", main = "Pop1")
lines(out_intervention[,3], type = "l", col = "blue")
legend("topright", legend = c("Control", "Intervention"), lwd = 2, col = c("orange3", "blue"))
plot(out[,6], type = "l", col = "orange3", main = "Pop2")
lines(out_intervention[,6], type = "l", col = "blue")


parameters["tbit"] <- 1e100
parameters["f12"] <- 0.003
parameters["f21"] <- 0.001  
parameters["beta2"] <- 0.25
out_high_migr <- deSolve::ode(y = istate, times = tps, func = SIR_patch, parms = parameters) 
parameters["tbit"] <- 10
out_intervention_high_migr <- deSolve::ode(y = istate, times = tps, func = SIR_patch, parms = parameters) 
par(mfrow=c(2,1))
plot(out_high_migr[,3], type = "l", col = "orange3", main = "Pop1")
lines(out_intervention_high_migr[,3], type = "l", col = "blue")
legend("topright", legend = c("Control", "Intervention"), lwd = 2, col = c("orange3", "blue"))
plot(out_high_migr[,6], type = "l", col = "orange3", main = "Pop2")
lines(out_intervention_high_migr[,6], type = "l", col = "blue")


