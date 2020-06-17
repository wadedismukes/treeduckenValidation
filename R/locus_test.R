

compare_locustree_mean_leaf_number <- function(gene_birth_rate,
                                     gene_death_rate,
                                     fixed_depth){
    speciesTree <- treeducken::sim_sptree_bdp_time(sbr = 1,
                                                   sdr = 0.25,
                                                   numbsim = 1,
                                                   t = fixed_depth)
    numberSpeciesTips <- length(speciesTree$tip.label)
    locusTreeSet <- treeducken::sim_locustree_bdp(species_tree = speciesTree[[1]],
                                  gbr = gene_birth_rate,
                                  gdr = gene_death_rate,
                                  lgtr = 0.0,
                                  num_loci = 10000)
    expectedValue <- treeducken::calculate_expected_leaves_locustree(t = fixed_depth,
                                                         dup_rate = gene_birth_rate,
                                                         loss_rate = gene_death_rate,
                                                         num_species = numberSpeciesTips)
    vectorOfNumberLoci <- vector(length = 10000)
    for(i in 1:10000){
        num_extant <- length(locusTreeSet[[i]]$tip.label) - length(grep("X", locusTreeSet[[i]]$tip.label))
        vectorOfNumberLoci[i] <- num_extant
    }
    print(mean(vectorOfNumberLoci) - expectedValue)
    (mean(vectorOfNumberLoci) - expectedValue) / sqrt(var(vectorOfNumberLoci))
}
