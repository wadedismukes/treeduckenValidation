gene_tree_msc_test <- function(sbr, gt_reps, test_reps) {
    # draw uniform (20,200)
    tips <- runif(test_reps, 20, 200)
    # sim sptree with GSA set sbr
    sptree <- treeducken::sim_sptree_bdp(sbr = 1.0,
                                         sdr = 0.0,
                                         num_tips = tips,
                                         numbsim = test_reps)
    # draw theta from lognormal(14,0.4) # ???
    theta <- stats::rlnorm(n = test_reps, meanLog = 14, sdlog = 0.4)
    tmrcas <- matrix(nrow = gt_reps, ncol = test_reps)
    for(i in 1:test_reps) {
        gts <- treeducken::sim_multispecies_coal(species_tree = sptree[[i]],
                                                 ne = theta[i],
                                                 num_sampled_individuals = 1,
                                                 num_genes = gt_reps)
        gts_summary_stats <- treeducken::genetree_summary_stat(gts,
                                                         locus_tree_indx = 1)
        tmrcas[, i] <- gts_summary_stats$TMRCA
    }
    # sim `reps` genetrees
    # calculate genetree_summary stats based on those trees
    # null the trees (to avoid memory overload)
    # return dataframe with the summary stats
    # this should be a dataframe with 10000 rows and 8(?) columns
    # of different summary stats for genetrees
}

gene_tree_msc_plot <- function(msc_df, expectation) {
    # plot the mdc_df (TMRCA with the expectation
}

gene_tree_mlc_test <- function(sbr,
                               popsize,
                               gbr,
                               locus_tree_reps,
                               gene_tree_reps,
                               test_reps) {
    # draw uniform(20, 200)
    tips <- runif(test_reps, 20, 200)

    # sim sptree with GSA set sbr
    sptree <- treeducken::sim_sptree_bdp(sbr = 1.0,
                                         sdr = 0.0,
                                         num_tips = tips,
                                         numbsim = test_reps)
    for(i in seq_len(length(sptree))) {
    # simulate locus tree with gbr locus_tree_reps times
        loc_trees <- treeducken::sim_locustree_bdp(species_tree = sptree[[i]],
                                                gbr = gbr,
                                                gdr = 0.0,
                                                lgtr = 0.0,
                                                num_loci = locus_tree_reps)
        # split each loocus tree into subtrees
        
    }
    # split each loocus tree into subtrees
    # simulate msc on each subtree with popsize set
    # calculate summary stats for each locus tree set
    # of subtrees compare with the expectation
}

