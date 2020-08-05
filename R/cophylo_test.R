host_tree_test <- function(hbr, hdr, cosp_rate, reps, time) {
    # simulate reps with treeducken
    rates_vec <- cbind(hbr, hdr)
    host_tree_num_tips <- matrix(nrow = 10000, ncol = rates_vec)
    for(i in seq_len(nrow(rates_vec))) {
        host_tree_pairs <- treeducken::sim_cophylo_bdp(hbr = rates_vec[i, 1],
                                                   hdr = rates_vec[i, 2],
                                                   cosp_rate = cosp_rate,
                                                   sbr = 1.0,
                                                   sdr = 0.5,
                                                   host_exp_rate = 0.0,
                                                   numbsim = reps,
                                                   time_to_sim = time)
        h_trees <- treeducken::host_tree(host_tree_pairs)
        # calculate number of host tips
        host_tree_num_tips[, i] <- ape::Ntip.multiPhylo(h_trees)
    }
    # return a dataframe of host tips (reps x (length(hbr) + length(hdr))
    host_tree_num_tips
}

host_tree_plot <- function(host_tip_df, expected_tips) {
    # plot of different reps with the lines for the expected
}

# for symb tree hbr = 0, hdr = 0
symb_tree_test <- function(cosp_rate, sbr, sdr, reps, time) {
    # simulates reps with treeducken
    # calculate number of tips
    rates_vec <- cbind(sbr, sdr)
    symb_tree_num_tips <- matrix(nrow = 10000, ncol = rates_vec)
    for(i in seq_len(nrow(rates_vec))) {
        tree_pairs <- treeducken::sim_cophylo_bdp(hbr = 0.0,
                                                   hdr = 0.0,
                                                   cosp_rate = 1.0,
                                                   sbr = rates_vec[i, 1],
                                                   sdr = rates_vec[i, 2],
                                                   host_exp_rate = 0.0,
                                                   numbsim = reps,
                                                   time_to_sim = time)
        s_trees <- treeducken::symb_tree(tree_pairs)
        # calculate number of tips
        symb_tree_num_tips[, i] <- ape::Ntip.multiPhylo(s_trees)
    }
    # return a dataframe of symb tips (reps x (length(sbr) + length(sdr))
    symb_tree_num_tips
}

symb_tree_plot <- function(symb_tip_df, expected_tips) {
    # same plot as host tree basically
}