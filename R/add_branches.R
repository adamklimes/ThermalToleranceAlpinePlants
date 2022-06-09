# add_branches
library(ape)
library(phytools)

# function which creates random bifurcations of selected tip-branch
# phylo: phylogenetic tree
# n: number of new tips
# tip_name: tip.label of the tip which should bifurcate
# new_names: optional vector of tip.labels for new tips

add_branches <- function(phylo, n = 1, tip_name, new_names = NULL){
  if (missing(new_names)) new_names <- paste0("new_t", 1:n) 
  pos <- which(phylo$tip.label == tip_name)
  max_len <- phylo$edge.length[which(phylo$edge[,2] == pos)]
  vec_dist <- sort(runif(n, 0, max_len), decreasing = TRUE)
  sel_tip_name <- tip_name
  for (i in 1:n){
    pos <- which(phylo$tip.label == sel_tip_name)
    phylo <- bind.tip(phylo, new_names[i], edge.length = vec_dist[i], where = pos, position = vec_dist[i])
    sel_tip_name <- sample(c(tip_name, new_names[1:i]), 1)
  }
  phylo
}

#tree1 <- pbtree(n = 10, tip.label = LETTERS[1:10])
#tree2 <- add_branches(tree1, 5, "A")
#plot(tree1)
#plot(tree2)