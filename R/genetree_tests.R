compare_genetree_tmrca_msc <- function(num_tips_species_tree, reps, theta){
    for(i in 1:reps){
        speciesTree <- treeducken::sim_sptree_bdp(n_tips = num_tips_species_tree,
                                                  sbr = 1,
                                                  sdr = 0.25,
                                                  numbsim = 1)
        geneBirthRate <- 0
        geneDeathRate <- 0
        transferRate <- 0
        numLoci <- 1
        numSampledIndividuals <- 1
        numGenesPerLocus <- reps
        locus_genetree <- treeducken::sim_locustree_genetree_mlc(species_tree = speciesTree[[1]],
                                                                 gbr = geneBirthRate,
                                                                 gdr = geneDeathRate,
                                                                 lgtr = transferRate,
                                                                 num_sampled_individuals = numSampledIndividuals,
                                                                 num_loci = 1,
                                                                 theta = theta,
                                                                 num_genes_per_locus = reps)
        stats_df <- treeducken::genetree_summary_stat(locus_tree_gene_tree_obj = locus_genetree,
                                                      locus_tree_indx = 1)
    }
    speciesTree <- treeducken::sim_sptree_bdp(n_tips = num_tips_species_tree,
                                              sbr = 1,
                                              sdr = 0.25,
                                              numbsim = 1)
    geneBirthRate <- 0
    geneDeathRate <- 0
    transferRate <- 0
    numLoci <- 1
    numSampledIndividuals <- 1
    numGenesPerLocus <- reps
    locus_genetree <- treeducken::sim_locustree_genetree_mlc(species_tree = speciesTree[[1]],
                                                             gbr = geneBirthRate,
                                                             gdr = geneDeathRate,
                                                             lgtr = transferRate,
                                                             num_sampled_individuals = numSampledIndividuals,
                                                             num_loci = 1,
                                                             theta = theta,
                                                             num_genes_per_locus = reps)
    stats_df <- treeducken::genetree_summary_stat(locus_tree_gene_tree_obj = locus_genetree,
                                                  locus_tree_indx = 1)
    stats_df
}
