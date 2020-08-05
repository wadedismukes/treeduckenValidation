library(treeducken)
library(treeduckenValidation)

reps <- 4
speciation_rates <- c(2.0, 4.0)
extinction_rates <- speciation_rates / 2
sim_age <- c(2.0, 4.0)
speciation_test_2_df <- matrix(nrow = reps, ncol = 8)
k <- 1
for(i in 1:length(speciation_rates)){
    for(j in 1:length(sim_age)) {
        speciation_test_2_df[1:reps,k] <- get_sptree_leaf_number_dist_time(speciation_rates[i],
                                                                extinction_rates[i],
                                                                time = sim_age[j],
                                                                reps = reps)
        k <- k + 1
        treesim_trees <- TreeSim::sim.bd.age(age = sim_age[j],
                                              lambda = speciation_rates[i],
                                              mu = extinction_rates[i],
                                              numbsim = reps) 
        for(t in 1:length(treesim_trees)) {
            speciation_test_2_df[t, k] <- length(treesim_trees[[t]][2])
        }
        k <- k + 1
    }
}
treesim_trees <- NULL
