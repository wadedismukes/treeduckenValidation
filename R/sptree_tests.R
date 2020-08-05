gsa_test_data_generation <- function(sbr, sdr, number_tips, reps) {
    # simulate reps with treeducken
    rates_vec <- cbind(sbr, sdr)
    gsa_sptree_tips <- matrix(nrow = 10000, ncol = rates_vec * 2)
    for(i in seq_len(nrow(rates_vec))) {
        td_species_trees <- treeducken::sim_sptree_bdp(sbr = rates_vec[i, 1],
                                                   sdr = rates_vec[i, 2],
                                                   n_tips = number_tips,
                                                   numbsim = reps,
                                                   time_to_sim = time)
        ts_species_trees <- TreeSim::sim.bd.taxa(n = number_tips,
                                                numbsim = reps,
                                                lambda =  rates_vec[i, 1],
                                                mu = rates_vec[i, 2])
        # calculate number of host tips
        gsa_sptree_tips[, i] <- ape::Ntip.multiPhylo(td_species_trees)
        gsa_sptree_tips[, i + nrow(rates_vec)] <- ape::Ntip.multiPhylo(ts_species_trees)
    }
    # return a dataframe of host tips (reps x (length(hbr) + length(hdr))
    gsa_sptree_tips
}

gsa_test_plot <- function(gsa_test_dataframe) {
    # plot boxplots? with the TreeSim plots next to the treeducken ones
}


ssa_test_data_generation <- function(sbr, sdr, time, reps) {
    # simulate the reps for the 10 pairs of birth and death rates
    # calculate number of tips
    # return a dataframe (10 columns, 10000 rows)
    rates_vec <- cbind(sbr, sdr)
    ssa_sptree_node_depths <- matrix(nrow = 10000, ncol = rates_vec)
    for(i in seq_len(nrow(rates_vec))) {
        trees <- treeducken::sim_sptree_bdp_time(sbr = rates_vec[i, 1],
                                                 sdr = rates_vec[i, 2],
                                                 numbsim = reps,
                                                 t = time)
        # calculate number of host tips
        ssa_sptree_node_depths[, i] <- max(ape::node.depth.edgelength(trees))
    }
    # return a dataframe of host tips (reps x (length(hbr) + length(hdr))
    ssa_sptree_node_depths
}

ssa_test_plot <- function(ssa_test_dataframe, expectation) {
    # plot with lines indicating the expectation
}