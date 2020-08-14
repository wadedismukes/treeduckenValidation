locus_tree_test <- function(gbr, gdr, number_tips, reps) {
    # simulate reps for 25 different combos of gbr and gdr
    # these are input as 2 vectors each of length 5
    # sbr and sdr don't change, time doesnt change
# for symb tree hbr = 0, hdr = 0
    # simulates reps with treeducken
    # calculate number of tips
    rates_vec <- cbind(gbr, gdr)
    loc_tree_num_tips <- matrix(nrow = reps + 1, ncol = nrow(rates_vec))
    for(i in seq_len(nrow(rates_vec))) {
        sp_tree <- treeducken::sim_sptree_bdp(sbr = 1.0,
                                              sdr = 0.0,
                                              numbsim = 1,
                                              n_tips = number_tips)
        loc_trees <- treeducken::sim_locustree_bdp(species_tree = sp_tree[[1]],
                                                    gbr = rates_vec[i, 1],
                                                    gdr = rates_vec[i, 2],
                                                    lgtr = 0.0,
                                                    num_loci = reps)

        # calculate number of tips
        loc_tree_num_tips[1:reps, i] <- ape::Ntip.multiPhylo(loc_trees)
        loc_tree_num_tips[reps + 1, i] <- max(ape::node.depth.edgelength(sp_tree[[1]]))
    }
    # return a dataframe of symb tips (reps x (length(sbr) + length(sdr))
    loc_tree_num_tips
}
    # returns a dataframe with 25 columns and 10000 rows


format_loctree_test <- function(loctr_tips, gbr) {
    #
    loctr_tips_df <- as.data.frame(loctr_tips)
    loc_column_namer <- function(gbr) {
        paste("delta=", gbr)
    }
    colnames(loctr_tips_df) <- loc_column_namer(gbr)
    loctr_tips_tidy_df <- tidyr::gather(loctr_tips_df)

}

get_loctr_expected_values <- function(gbr,
                                      gdr,
                                      num_tips,
                                      loctr_tips,
                                      loctr_tips_tidy_df) {
    loct_expected_tips <- vector(length = length(gbr))
    for(i in seq_len(length(gbr))) {
        loct_expected_tips[i] <- treeducken::calculate_expected_leaves_locustree(
            t = loctr_tips[reps + 1, i],
            dup_rate = gbr[i],
            loss_rate = gdr[i],
            num_species = num_tips)
    }

    loct_expected_tips_df <- as.data.frame(t(loct_expected_tips))
    colnames(loct_expected_tips_df) <- unique(loctr_tips_tidy_df[,1])
    tidyr::gather(loct_expected_tips_df)
}

