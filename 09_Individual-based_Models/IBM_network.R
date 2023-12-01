library("igraph")
init_S <- 100
init_I <- 10
beta <- 0.5
gamma <- 1 / 5
tau <- 1 / 10
alpha <- 1
no_of_timesteps <- 100
pop <- c(rep(0, init_S), rep(1, init_I))
sim_table <- as.data.frame(matrix(NA, no_of_timesteps, 5))
colnames(sim_table) <- c("Day", "S", "E", "I", "R")


connectivity <- 2
rewiring <- 0.15
nw <- sample_smallworld(1, length(pop), connectivity, rewiring)
#mean_distance(nw)
#transitivity(nw, type="average")
plot(nw)
as_edgelist(nw, names = TRUE)
neighborhood(nw, nodes = 1:3)


for(t in 1:no_of_timesteps){ 

	## Simulate what happens at each timestep
	sim_table$Day[t] <- t
	sim_table$S[t] <- sum(pop == 0)
	sim_table$E[t] <- sum(pop == 1)
	sim_table$I[t] <- sum(pop == 2)
	sim_table$R[t] <- sum(pop == 3)

	
	random_gamma <- runif(length(pop))
	random_tau <- runif(length(pop))
	random_alpha <- runif(length(pop))
  
	for(i in seq_along(pop)){
    	#Simulate what happens to each individual
    
		## INFECTION ##
		if(pop[i] == 2) {
		# how many infectious contacts do you have today
		ncontacts <- rpois(1, beta)
			if (ncontacts > 0) {
				# who do you have a contact with?
				pcontacted <- sample(as.vector(neighborhood(nw, nodes = i)[[1]][-1]), ncontacts, replace = T)
				#print(paste("p",ncontacts))
				#print(paste("w",pcontacted))
				for (w in 1:ncontacts){
					if (pop[pcontacted[w]] == 0) {
						pop[pcontacted[w]] <-1
					}
				}   
			}
		} 
	
		## LATENCY ##
		if(pop[i] == 1){
			if(random_gamma[i] < gamma) { 
				pop[i] <- 2
			}
		}

		## RECOVERY ##
		if(pop[i] == 2){
			if(random_tau[i] < tau) {
				pop[i] <- 3
			}
		}
	
		## LOSS OF IMMUNITY ##
		if(pop[i] == 3){
			if(random_alpha[i] < alpha) {
				pop[i] <- 0
			}
		}
    
	} # close loop for individuals

} # close loop for timesteps

# look at the results table
sim_table

# plot results
plot(sim_table$Day, sim_table$S, type = "l", ylim = c(0, length(pop)),
    col = "black", xlab = "Weeks", ylab = ("# People"))
lines(sim_table$Day, sim_table$I, col = "red")
lines(sim_table$Day, sim_table$R, col = "blue")
legend("topright", legend = c("S", "I", "R"), lwd = 1,
    col = c("black", "red", "blue"))
