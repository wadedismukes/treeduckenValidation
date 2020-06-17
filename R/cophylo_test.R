# host tree
#
#
compare_host_tree_with_expectation <- function(hbr, hdr, cosp_rate,
                                       sbr, sdr, host_exp_rate,
                                       num_replicates = 10000, time_to_sim){
    host_tree_pairs <- treeducken::sim_cophylo_bdp(hbr = hbr,
                                                   hdr = hdr,
                                                   cosp_rate = cosp_rate,
                                                   sbr = sbr,
                                                   sdr = sdr,
                                                   host_exp_rate = host_exp_rate,
                                                   numbsim = num_replicates,
                                                   time_to_sim = time_to_sim)
    expectedValue <- treeducken::calculate_expected_leaves_sptree(lambda = hbr,
                                                                  mu = sdr,
                                                                  t = time)
    vectorOfNumberSpecies <- vector(length = 10000)
    for(i in 1:10000){
        num_extant <- length(host_tree_pairs[[i]]$host_tree$tip.label)
        vectorOfNumberSpecies[i] <- num_extant
    }
    (mean(vectorOfNumberSpecies) - expectedValue) / sqrt(var(vectorOfNumberSpecies))
}
