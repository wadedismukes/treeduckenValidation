gsa_test_data_generation <- function(sbr, sdr, number_tips, reps) {
    # simulate reps with treeducken
    rates_vec <- cbind(sbr, sdr)
    gsa_sptree_node_depths <- matrix(nrow = reps, ncol = nrow(rates_vec) * 2)
    for(i in seq_len(nrow(rates_vec))) {
        td_species_trees <- treeducken::sim_sptree_bdp(sbr = rates_vec[i, 1],
                                                   sdr = rates_vec[i, 2],
                                                   n_tips = number_tips,
                                                   numbsim = reps)
        ts_species_trees <- TreeSim::sim.bd.taxa(n = number_tips,
                                                numbsim = reps,
                                                lambda =  rates_vec[i, 1],
                                                mu = rates_vec[i, 2])
        # calculate number of host tips
        multiphy_calc_edges <- function(x) max(ape::node.depth.edgelength(x))
        gsa_sptree_node_depths[, i] <- sapply(td_species_trees, multiphy_calc_edges)
        gsa_sptree_node_depths[, i + nrow(rates_vec)] <- sapply(ts_species_trees, multiphy_calc_edges)

    }
    # return a dataframe of host tips (reps x (length(hbr) + length(hdr))
    gsa_sptree_node_depths
}

format_gsa_test_df <- function(gsa_spt_nds, sbr) {
    # plot boxplots? with the TreeSim plots next to the treeducken ones
    gsa_spt_nds_df <- as.data.frame(gsa_spt_nds)
    column_namer <- function(spec_rate, sim = "Treeducken") {
        paste("lambda=", spec_rate, ",", sim, sep="")
    }
    colnames(gsa_spt_nds_df) <- c(column_namer(sbr),
                                  column_namer(sbr, sim = "TreeSim"))
    tidy_df <- tidyr::gather(gsa_spt_nds_df)
    split_names <- stringr::str_split(tidy_df$key, ",", simplify = TRUE)
    tidy_df$key <- split_names[,1]
    tidy_df$simulator <- split_names[,2]
    tidy_df
}


ssa_test_data_generation <- function(sbr, sdr, time, reps) {
    # simulate the reps for the 10 pairs of birth and death rates
    # calculate number of tips
    # return a dataframe (10 columns, 10000 rows)
    rates_vec <- cbind(sbr, sdr)
    ssa_sptree_node_depths <- matrix(nrow = reps, ncol = nrow(rates_vec))
    for(i in seq_len(nrow(rates_vec))) {
        trees <- treeducken::sim_sptree_bdp_time(sbr = rates_vec[i, 1],
                                                 sdr = rates_vec[i, 2],
                                                 numbsim = reps,
                                                 t = time)
        # calculate number of host tips
        ssa_sptree_node_depths[, i] <- ape::Ntip.multiPhylo(trees)

    }
    # return a dataframe of host tips (reps x (length(hbr) + length(hdr))
    ssa_sptree_node_depths
}

format_ssa_test <- function(ssa_spt_tips, sbr) {
    ssa_spt_tips_df <- as.data.frame(ssa_spt_tips)
    ssa_column_namer <- function(spec_rate) {
        paste("lambda=", spec_rate)
    }
    colnames(ssa_spt_tips_df) <- ssa_column_namer(sbr)
    tidyr::gather(ssa_spt_tips_df)
}

get_ssa_expected_values <- function(sbr, sdr, sim_time, ssa_spt_tips_tidy_df) {
    ssa_expected_tips <- vector(length = length(sbr))
    for(i in seq_len(length(sbr))) {
        ssa_expected_tips[i] <- treeducken::calculate_expected_leaves_sptree(sbr[i],
                                                                             sdr[i],
                                                                             sim_time)
    }
    ssa_exp_tips_df <- as.data.frame(t(ssa_expected_tips))
    colnames(ssa_exp_tips_df) <- unique(ssa_spt_tips_tidy_df[,1])
    tidyr::gather(ssa_exp_tips_df)

}