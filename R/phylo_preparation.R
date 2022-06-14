# phylo_preparation
library(ape)
library(Taxonstand)
library(V.PhyloMaker)

# data_loading
dat <- read.csv("data/Dataset_S1.csv")

# data_preparation
# TPL(unique(dat$Species), file = "data/analyses/taxonstandNames.csv")
dat_taxonstand <- read.csv("data/analyses/taxonstandNames.csv")
  # use the original name for "Cerastium arvensiforme"
dat_taxonstand$New.Species[dat_taxonstand$Taxon == "Cerastium arvensiforme"] <- "arvensiforme"
# write.csv(dat_taxonstand, file = "data/analyses/taxonstandNames.csv", row.names = FALSE)
spec_clad <- with(dat_taxonstand, 
  data.frame(species = paste(New.Genus, New.Species),
    genus = New.Genus, family = Family))
spec_clad$family[spec_clad$species == "Cerastium sp."] <- "Caryophyllaceae"
spec_clad$family[spec_clad$species == "Silene rimbachii"] <- "Caryophyllaceae"
spec_clad$family[spec_clad$family == "Compositae"] <- "Asteraceae"
set.seed(61)
phy <- phylo.maker(spec_clad, scenarios = "S2")$scenario.2$run.1
# write.tree(phy, "data/phy.tre")

#_