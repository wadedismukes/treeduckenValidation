gene_tree_msc_test <- function(theta, tips, gt_reps, test_reps) {

    # sim sptree with GSA set sbr
    for(i in seq_len(length(tips))) {
        sptree <- treeducken::sim_sptree_bdp(sbr = 1.0,
                                             sdr = 0.0,
                                             numbsim = test_reps,
                                             n_tips = tips[i])
    }

    # draw theta from lognormal(14,0.4) # ???
    tmrcas <- matrix(nrow = gt_reps, ncol = test_reps)
    for(i in 1:test_reps) {
        gts <- treeducken::sim_multispecies_coal(species_tree = sptree[[i]],
                                                 ne = theta[i],
                                                 num_sampled_individuals = 1,
                                                 mutation_rate = 1e-9,
                                                 generation_time = 1e-6,
                                                 num_genes = gt_reps,
                                                 rescale = TRUE)
        gts_summary_stats <- treeducken::genetree_summary_stat(gts,
                                                         container_tree_indx = 1)
        tmrcas[, i] <- gts_summary_stats$tmrca
    }
    tmrcas
    # sim `reps` genetrees
    # calculate genetree_summary stats based on those trees
    # null the trees (to avoid memory overload)
    # return dataframe with the summary stats
    # this should be a dataframe with 10000 rows and 8(?) columns
    # of different summary stats for genetrees
}

format_gene_tree_df <- function(msc_tmrcas, popsize) {
    msc_gt_tmrca_df <- as.data.frame(msc_tmrcas)
    msc_column_namer <- function(theta) {
        paste("popsize=", theta)
    }
    colnames(msc_gt_tmrca_df) <- msc_column_namer(popsize)
    tidyr::gather(msc_gt_tmrca_df)
}

gene_tree_mlc_test <- function(sbr,
                               popsize,
                               gbr,
                               locus_tree_reps,
                               gene_tree_reps,
                               test_reps) {
    # draw uniform(20, 200)
    tips <- 50
    # sim sptree with GSA set sbr
    # sim sptree with GSA set sbr
    sptree <- treeducken::sim_sptree_bdp(sbr = 1.0,
                                        sdr = 0.0,
                                        n_tips = tips,
                                        numbsim = test_reps)
    tmrcas <- matrix(nrow = gene_tree_reps, ncol = locus_tree_reps * test_reps)
    for(i in seq_len(test_reps)) {
        # simulate locus tree with gbr locus_tree_reps times
        loc_trees <- treeducken::sim_locustree_bdp(species_tree = sptree[[i]],
                                                gbr = gbr,
                                                gdr = 0.0,
                                                lgtr = 0.0,
                                                num_loci = locus_tree_reps)
        # split each loocus tree into subtrees

        for(j in seq_len(locus_tree_reps)) {
            gene_trees <- treeducken::sim_multilocus_coal(
                locus_tree = loc_trees[[j]],
                effective_pop_size = popsize[i],
                mutation_rate = 1e-9,
                generation_time = 1e-6,
                num_reps = gene_tree_reps)
            mtrees <- treeducken::retrieve_parent_genetrees(gene_trees)
            col_indx <- j + (test_reps * (i - 1))
            tmrcas[, col_indx] <- unlist(
                                        lapply(mtrees,
                                            function(x)
                                            max(ape::node.depth.edgelength(x))))
        }
    }
    tmrcas
}

