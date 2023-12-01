setwd("C:/mgh/IBM")
library("plot.matrix")
library("igraph")

init_S <- 5000
init_I <- 500
p <- 0.1
gamma <- 1 / 5
tau <- 1 / 10
alpha <- 1
no_of_timesteps <- 100

pop <- rep(0, init_S + init_I)
sim_table <- as.data.frame(matrix(NA, no_of_timesteps, 5))
colnames(sim_table) <- c("Day", "S", "E", "I", "R")

## import contact matrices
c_home <- as.matrix(read.table("contact_home.txt", sep = ",", header = F))
c_school <- as.matrix(read.table("contact_school.txt", sep = ",", header = F))
c_work <- as.matrix(read.table("contact_work.txt", sep = ",", header = F))
c_other <- as.matrix(read.table("contact_other.txt", sep = ",", header = F))
total_contacts <- c_home + c_school + c_other
plot(total_contacts)

## import age structure
popstruc <- read.csv("pop_structure_X.csv", header = T)
# convert from 1000s to total numbers
popstruc[, 2] <- popstruc[, 2] * 1000

## assign an age group to each individual
samplesize <- length(pop)
ages <- with(popstruc, sample(length(popstruc[, 1]),
    samplesize, p = (popstruc[, 2] / sum(popstruc[, 2])), replace = T)) 

## build network based on age
X <- as.matrix(total_contacts)
g <- graph_from_adjacency_matrix(X, mode = "lower", 
    weighted = "weight", diag = T)
g2 <- delete.edges(g, which(E(g)$weight < 1))
plot(g2, edge.width = E(g2)$weight)

## init_Ialise infections in random individuals
starting_infs <- sample(length(pop), init_I)
pop[starting_infs] <- 1

for (t in 1:no_of_timesteps) {
    #Simulate what happens at each timestep
    sim_table$Day[t] <- t
    sim_table$S[t] <- sum(pop == 0)
    sim_table$E[t] <- sum(pop == 1)
    sim_table$I[t] <- sum(pop == 2)
    sim_table$R[t] <- sum(pop == 3)

    random_gamma <- runif(length(pop))
    random_tau <- runif(length(pop))
    random_alpha <- runif(length(pop))

    for (i in seq_along(pop)) {
    #Simulate what happens to each individual

        ## INFECTION ##
        if (pop[i] == 2) {
            # how many infectious contacts do you have today
            ncontacts <- rpois(1, p * sum(total_contacts[ages[i], ]))
            # print(ncontacts)
            if (ncontacts > 0) { 
                for (w in 1:ncontacts) {
                    # which age group do you get in contact with
                    target_age <- sample(length(total_contacts[, 1]), 1,
                    prob = total_contacts[ages[i], ], replace = TRUE)
                    target_person <- sample(which(ages == target_age), 1)
                    if (pop[target_person] == 0) {
                        pop[target_person] <- 1
                    }
                }
            }
        }

        ## LATENCY ##
        if (pop[i] == 1) {
            if (random_gamma[i] < gamma){ 
                pop[i] <- 2
            }
        }

        ## RECOVERY ##
        if(pop[i] == 2) {
            if(random_tau[i] < tau){
                pop[i] <- 3
            }
        }

        ## LOSS OF IMMUNITY ##
        if(pop[i] == 3) {
            if(random_alpha[i] < alpha) {
                pop[i] <- 0
            }
        }
    } # close loop for individuals
} # close loop for timesteps

# look at the results table
sim_table

plot(sim_table$Day, sim_table$S, type = "l", ylim = c(0, length(pop)),
    col = "black", xlab = "Weeks", ylab = ("# People"))
lines(sim_table$Day, sim_table$I, col = "red")
lines(sim_table$Day, sim_table$R, col = "blue")
legend("topright", legend = c("S", "I", "R"), lwd = 1,
    col = c("black", "red", "blue"))
