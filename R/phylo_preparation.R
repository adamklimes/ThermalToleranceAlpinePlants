# phylo_preparation
library(ape)

# data loading
dat <- read.table("data/Temperature_resistance_data.csv", sep = ";", header = TRUE)
spec_gform <- read.table("data/Growth_Forms_List.csv", sep = ";", header = TRUE)
ztree <- read.tree("data/Vascular_Plants_rooted.dated.tre") # phylogenetic tree from Zanne et al. 2014

# data preparation
dat$rng <- dat$He_lt50 - dat$Fr_lt50
dat$Species[dat$Species == "Pedicularis_parryi_late_melt"] <- "Pedicularis_parryi2"

# creating random tree
pom_spec <- unique(unlist(lapply(strsplit(dat$Species, "_"), function(x) paste(x[1], x[2], sep="_"))))
pom_spec[grepl("\"", pom_spec)] <- paste0(pom_spec[grepl("\"", pom_spec)], 2)
pom_spec <- gsub("\"", "", pom_spec)
old_spec <- c("Loricaria_antisanensis", "Oritrophium_peruvianum", "Oritrophium_hieracioides", "Lasiocephalus_ovatus", "Pentacalia_andicola")
new_spec <- c("Gamochaeta_sp", "Diplostephium_sp", "Diplostephium_sp2", "Senecio_sp", "Senecio_sp2")
pom_spec_new <- pom_spec
pom_spec_new[pom_spec_new %in% old_spec] <- new_spec

ztree2 <- add_branches(ztree, 5, "Werneria_nubigena", paste0("Werneria_sp",1:5))
ztree3 <- add_branches(ztree2, 5, "Bartsia_alpina", paste0("Bartsia_alsp",1:5))
ztree4 <- add_branches(ztree3, 5, "Bartsia_trixago", paste0("Bartsia_trsp",1:5))
ztree5 <- add_branches(ztree4, 5, "Diplostephium_rupestre", paste0("Diplostephium_newsp",1:5))
r_tree <- create_subtree(pom_spec_new, ztree5, randomSpec = TRUE, ntrees = 1)[[3]][[1]]
# save(r_tree, file="data/analyses/r_tree.RData")

#_