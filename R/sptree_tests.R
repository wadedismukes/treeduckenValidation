get_sptree_leaf_number_dist_time <- function(sbr, sdr, time, reps){
    speciesTrees <- treeducken::sim_sptree_bdp_time(sbr = sbr,
                                                   sdr = sdr,
                                                   numbsim = reps,
                                                   t = time)
    expectedValue <- treeducken::calculate_expected_leaves_sptree(lambda = sbr,
                                                                  mu = sdr,
                                                                  t = time)
    vectorOfNumberSpecies <- vector(length = 10000)
    for(i in 1:10000){
        num_extant <- length(speciesTrees[[i]]$tip.label)
        vectorOfNumberSpecies[i] <- num_extant
    }
    vectorOfNumberSpecies
}

get_sptree_leaf_number_dist <- function(sbr, sdr, num_tips, reps){
    speciesTrees <- treeducken::sim_sptree_bdp(sbr = sbr,
                                               sdr = sdr,
                                               numbsim = reps,
                                               n_tips = num_tips)
    expectedValue <- treeducken::estimate_node_heights(lambda = sbr,
                                                       mu = sdr,
                                                       n = num_tips)
    vectorOfNodeDepths <- vector(length = reps)
    for(i in 1:reps){
        nodeDepth <- max(ape::node.depth.edgelength(speciesTrees[[i]])) + speciesTrees[[i]]$root.edge
        vectorOfNodeDepths[i] <- nodeDepth
    }
    vectorOfNodeDepths
    #((mean(vectorOfNodeDepths) / sqrt(var(vectorOfNodeDepths))) - (mean(vectorOfNodeDepthsTreeSim)) / sqrt(var(vectorOfNodeDepthsTreeSim)))
    #(mean(vectorOfNodeDepths) - expectedValue) / sqrt(var(vectorOfNodeDepths))
}
