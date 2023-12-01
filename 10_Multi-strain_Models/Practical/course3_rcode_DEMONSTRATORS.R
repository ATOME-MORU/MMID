

require(MANTIS)##load the MANTIS R package

###########################################################
##BASE CODE FOR PART 1

#these next instructions define the time variables of the simulation
tMax= 700 #total time in years
tInt= 0.005 #simulation time step
tObsPer= 100 #output time step (it samples every tObsPer)

#these next instructions define the existing strains and their interactions
epiStruc= c(3,2) #antigenic structure: with 2 loci with 3 and 2 alleles each (6 strains)
nStrains= extractNumberStrains(epiStruc) #get number of strains in the system =3x2=6
infInitCond=c(0.000013, 0.000012, 0.00010, 0.000014, 0.000011, 0.000010) #initial infected of each strain
gamma= 0 #strength of cross-immunity between strains

#epidemiological parameters
beta= c(292,292,292,292,292,292) #transmission rates of the strains
sigma= (1/5)*365 #5 days infectious period
mu= 1/50 #50 years of host life-span

#this instruction executes the simulation
simdata= runMANTIS(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu)

#this next instruction plots the proportion of infected individuals
#notice the variables xiObs=0.6, xfObs=1.0, which define a "zoom in" window
#within 60 and 100% of the time range of the simulation
#you can change these numbers to zoom in at different times and explore the simulation

#the legend shows each strain, with its corresponding "sequence" of alleles
plotY(simdata, xiObs=0.6, xfObs=1.0, addLegend=TRUE)

###################################
## exercise 1.1

plotY(simdata, xiObs=0.05, xfObs=0.35, addLegend=TRUE)

###################################
## exercise 1.2b

extractYFinalConditions(simdata)




###################################
## exercise 1.3

beta= c(292,584,292,292,292,292) #transmission rates of the strains
simdata= runMANTIS(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu)
plotY(simdata, xiObs=0.05, xfObs=0.3, addLegend=TRUE)

###################################
## exercise 1.4

gamma<- 0
beta= c(292,292,292,292,292,292) #transmission rates of the strains
simdata= runMANTIS(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu)
plotYDiversity(simdata, xiObs=0.3, xfObs=1.0)

###################################
## exercise 2.1

gamma= 0.95 #strength of cross-immunity between strains
beta= c(292,292,292,292,292,292) #transmission rates of the strains
simdata= runMANTIS(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu)
plotY(simdata, xiObs=0.6, xfObs=1.0, addLegend=TRUE)

###################################
## exercise 2.2

extractYFinalConditions(simdata)

################################
## exercise 2.4

gamma= 0.95 #strength of cross-immunity between strains
beta= c(292,584,292,292,292,292) #transmission rates of the strains
simdata= runMANTIS(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu)
plotY(simdata, xiObs=0.6, xfObs=1.0, addLegend=TRUE)
extractYFinalConditions(simdata)


################################
## exercise 2.5

gamma= 0.95 #strength of cross-immunity between strains
beta= c(292,292,292,292,292,292) #transmission rates of the strains
simdata= runMANTIS(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu)
plotYDiversity(simdata, xiObs=0.3, xfObs=1.0)


################################
## exercise 3.1

gamma= 0.7 #strength of cross-immunity between strains
beta= c(292,292,292,292,292,292) #transmission rates of the strains
simdata= runMANTIS(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu)
plotY(simdata, xiObs=0.6, xfObs=0.65, addLegend=TRUE)
plotZ(simdata, xiObs=0.6, xfObs=0.65)
plotW(simdata, xiObs=0.6, xfObs=0.65)

################################
## exercise 3.4

newPlotW(simdata)
gamma= 0.7 #strength of cross-immunity between strains
beta= c(292,292,292,292,292,292) #transmission rates of the strains
mu = 1/70
simdata_human= runMANTIS(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu)
mu = 1/10
simdata_duck= runMANTIS(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu)


plotY(simdata_human, xiObs=0.6, xfObs=1.0, addLegend=TRUE)
plotY(simdata_duck, xiObs=0.6, xfObs=1.0, addLegend=TRUE)

################################
## exercise 3.5

gamma<- 0.7
mu = 1/70
beta= c(292,292,292,292,292,292) #transmission rates of the strains
simdata= runMANTIS(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu)
plotYDiversity(simdata, xiObs=0.3, xfObs=1.0)



###########################################################
##BASE CODE FOR PART 4

## scenario i - no control

tMax= 700
tInt= 0.005
tObsPer= 100 
epiStruc= c(3,2) 
nStrains= extractNumberStrains(epiStruc)
infInitCond=c(0.000013, 0.000012, 0.00010, 0.000014, 0.000011, 0.000010) 
gamma= 0.95 
beta= c(292,292,292,292,292,292) 
sigma= (1/5)*365
mu= 1/50 
simdata= runMANTIS(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu)
plotY(simdata, xiObs=0.6, xfObs=1.0, addLegend=TRUE)
extractYFinalConditions(simdata)

## scenario ii - control strain 2 reducing its beta to zero at time 500

tMax= 1000
tInt= 0.005
tObsPer= 100 
epiStruc= c(3,2) 
nStrains= extractNumberStrains(epiStruc)
infInitCond=c(0.000013, 0.000012, 0.00010, 0.000014, 0.000011, 0.000010) 
gamma= 0.95 
beta= c(292,292,292,292,292,292) 
sigma= (1/5)*365
mu= 1/50 
tBetaChange= 500 #time at which beta will be changed
changedStrain= 2 #which strain to change beta for
changedBeta= 292*0 #set the new beta
simdata<- runInvasionWithOneBetaChange(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu,
                                       tBetaChange, changedBeta, changedStrain) #run MANTIS
plotY(simdata, xiObs=0.45, xfObs=0.99, ymax=0.001, addLegend=TRUE) #plot entire simulation
extractYFinalConditions(simdata)

## scenario iii - control strain 6 reducing its beta to zero at time 500

tMax= 700
tInt= 0.005
tObsPer= 100 
epiStruc= c(3,2) 
nStrains= extractNumberStrains(epiStruc)
infInitCond=c(0.000013, 0.000012, 0.00010, 0.000014, 0.000011, 0.000010) 
gamma= 0.95 
beta= c(292,292,292,292,292,292) 
sigma= (1/5)*365
mu= 1/50
tBetaChange= 500 #time at which beta will be changed
changedStrain= 6 #which strain to change beta for
changedBeta= 292*0 #set the new beta
simdata<- runInvasionWithOneBetaChange(epiStruc, tMax, tObsPer, tInt, infInitCond, beta, gamma, sigma, mu,
                                       tBetaChange, changedBeta, changedStrain) #run MANTIS
plotY(simdata, xiObs=0.5, xfObs=0.99, ymax=0.001, addLegend=TRUE) #plot entire simulation
extractYFinalConditions(simdata)



