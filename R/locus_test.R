locus_tree_test <- function(gbr, gdr, time, reps) {
    # simulate reps for 25 different combos of gbr and gdr
    # these are input as 2 vectors each of length 5
    # sbr and sdr don't change, time doesnt change
# for symb tree hbr = 0, hdr = 0
    # simulates reps with treeducken
    # calculate number of tips
    rates_vec <- cbind(gbr, gdr)
    loc_tree_num_tips <- matrix(nrow = 10001, ncol = rates_vec)
    for(i in seq_len(nrow(rates_vec))) {
        sp_tree <- treeducken::sim_sptree_bdp_time(sbr = 1.0,
                                                   sdr = 0.5,
                                                   t = time,
                                                   numbsim = 1)
        loc_trees <- treeducken::sim_locustree_bdp(species_tree = sp_tree,
                                                    gbr = rates_vec[i, 1],
                                                    gdr = rates_vec[i, 2],
                                                    lgtr = 0.0,
                                                    num_loci = reps)
        # calculate number of tips
        loc_tree_num_tips[1:10000, i] <- ape::Ntip.multiPhylo(loc_trees)
        loc_tree_num_tips[10001, i] <- ape::Ntip.multiPhylo(sp_tree)
    }
    # return a dataframe of symb tips (reps x (length(sbr) + length(sdr))
    loc_tree_num_tips
}
    # returns a dataframe with 25 columns and 10000 rows


locus_tree_plot <- function(locustree_test_df, expected_value) {
    # plot boxplots of a dataframe with 
}



