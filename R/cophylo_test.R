host_tree_test <- function(hbr, hdr, cosp_rate, reps, time) {
    # simulate reps with treeducken
    rates_vec <- cbind(hbr, hdr)
    host_tree_num_tips <- matrix(nrow = reps, ncol = nrow(rates_vec))
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

format_host_test_data <- function(ht_test_tips, lambda_h) {
    ht_test_tips_df <- as.data.frame(ht_test_tips)
    ht_column_namer <- function(lambda_h) {
        paste("lambda_h=", lambda_h)
    }
    colnames(ht_test_tips_df) <- ht_column_namer(lambda_h)
    tidyr::gather(ht_test_tips_df)
}

get_host_expected_values <- function(lambda_h,
                                     lambda_c,
                                     mu_h,
                                     sim_time,
                                     ht_test_tidy_df) {
    ht_expected_tips <- vector(length = length(lambda_h))

    for(i in seq_len(length(lambda_h))) {
        ht_expected_tips[i] <- treeducken::calculate_expected_leaves_sptree(
            lambda = lambda_c + lambda_h[i],
            mu = mu_h[i],
            t = sim_time)
    }
    ht_expected_tips_df <- as.data.frame(t(ht_expected_tips))
    colnames(ht_expected_tips_df) <- unique(ht_test_tidy_df[,1])
    tidyr::gather(ht_expected_tips_df)
}

# for symb tree hbr = 0, hdr = 0
symb_tree_test <- function(cosp_rate, sbr, sdr, reps, time) {
    # simulates reps with treeducken
    # calculate number of tips
    rates_vec <- cbind(sbr, sdr)
    symb_tree_num_tips <- matrix(nrow = reps, ncol = nrow(rates_vec) + nrow(rates_vec))
    for(i in seq_len(nrow(rates_vec))) {
        tree_pairs <- treeducken::sim_cophylo_bdp(hbr = 0.0,
                                                   hdr = 0.0,
                                                   cosp_rate = cosp_rate,
                                                   sbr = rates_vec[i, 1],
                                                   sdr = rates_vec[i, 2],
                                                   host_exp_rate = 0.0,
                                                   numbsim = reps,
                                                   time_to_sim = time)
        h_trees <- treeducken::host_tree(tree_pairs)
        s_trees <- treeducken::symb_tree(tree_pairs)
        # calculate number of
        symb_tree_num_tips[, i + nrow(rates_vec)] <- ape::Ntip.multiPhylo(h_trees)
        symb_tree_num_tips[, i] <- ape::Ntip.multiPhylo(s_trees)
    }
    # return a dataframe of symb tips (reps x (length(sbr) + length(sdr))
    symb_tree_num_tips
}

symb_tree_plot <- function(symb_tip_df, expected_tips) {
    # same plot as host tree basically
}