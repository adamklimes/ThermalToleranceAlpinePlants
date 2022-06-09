# create_subtree function

# function which takes species and tree and output subtree of these species
# - it works also when species are missing in the tree:
#	function selects some other species from the tree from the same genus
#	if the genus is monophyletic in relation to other selected species
#	
#	In case of non-monophyletic genus function can pick species at random 
#	and create multiple such random subtrees.    

#	Also more than one species identification is possible to use 
#	(e.i.: the species is called X, if not found look for it under the name Y)

# species:  vector of names of target species; or list of several such vectors as additional identifications
# tree:	phylogeny from which subtree should be created
# assignGenus: T/F - if for missing species should be selected some other in case the genus is monophyletic 
# randomSpec: T/F - if for non-monophyletic random species from the genus should be selected
# ntrees: 	in case that randomSpec==T, how many such trees should be created

library(ape)
create_subtree <- function(species,tree,assignGenus=T,randomSpec=F,ntrees=100){
	sel.gen <- function(x) unlist(lapply(strsplit(x,"_"),function(y)y[1]))
	if (class(species)!="list") {species <- gsub(" ","_",species)
		prez.spec <- species %in% tree$tip
		spec <- species
		assign.id <- rep(NA,length(spec))
		assign.id[prez.spec] <- 1 
	} else {
		species <- lapply(species,function(x) gsub(" ","_",x))
		prez.spec <- rep(FALSE,length(species[[1]]))
		assign.id <- rep(NA,length(species[[1]]))
		for (i in 1:length(species)){
			n.prez.spec <- (species[[i]] %in% tree$tip)|(prez.spec)
			assign.id[which(prez.spec!=n.prez.spec)] <- i
			prez.spec <- n.prez.spec
		}
		spec <- species[[1]]		
	}
	genus.spec <- NULL	
	if (assignGenus){
		full.genus.names <- sel.gen(spec[!prez.spec])
		genus.names <- full.genus.names[!(full.genus.names %in% sel.gen(spec[prez.spec]))]
		genus.names <- genus.names[!(genus.names %in% genus.names[duplicated(genus.names)])]
		tree.genus <- sel.gen(tree$tip)
		genus.spec <- tree$tip[tree.genus %in% genus.names]
		genus.spec <- genus.spec[!(genus.spec %in% spec[prez.spec])]
		add.genus.spec <- tree$tip[tree.genus %in% full.genus.names]
		add.genus.spec <- add.genus.spec[!(add.genus.spec %in% genus.spec) & !(add.genus.spec %in% spec[prez.spec])]
	}
	if (class(species)!="list") {pom.spec <- spec[prez.spec]}else{
		pom.spec <- rep(NA,length(spec))
		for (i in 1:length(spec)){
			pom <- species[[assign.id[i]]][i]
			if (is.null(pom)) pom <- NA
			pom.spec[i] <- pom
		}
		pom.spec <- pom.spec[!is.na(pom.spec)]
		 
	}
	sel.spec<- c(pom.spec,genus.spec,add.genus.spec)
	tree.out <- drop.tip(tree,tree$tip[!(tree$tip %in% sel.spec)])
	non.mon.gen <- NULL
	if (assignGenus){
		tree.out.clear <- drop.tip(tree.out,add.genus.spec)
		d.mat <- cophenetic(tree.out.clear)
		d.mat <- d.mat[rownames(d.mat) %in% genus.spec,]
		colnames(d.mat) <- sel.gen(colnames(d.mat))
		found.genus <- unique(tree.genus[tree.genus %in% genus.names])
		d.mat <- t(apply(d.mat,1,rank))
		mon.gen <- rep(NA,nrow(d.mat))
		for (i in 1:nrow(d.mat)){
			sel.col <- colnames(d.mat)==(sel.gen(rownames(d.mat))[i]) 
			mon.gen[i] <- max(d.mat[i,sel.col])<=sum(sel.col)
		}
		pom <- sample(rownames(d.mat)[mon.gen])
		genus.to.dump <- pom[duplicated(sel.gen(pom))]
		non.mon.gen <- spec[sel.gen(spec) %in% sel.gen(rownames(d.mat)[!mon.gen])]
		if (randomSpec) {tree.out.add <- drop.tip(tree.out,genus.to.dump)
			add.genus.spec <- c(add.genus.spec, rownames(d.mat)[!mon.gen])
		}
		tree.out <- drop.tip(tree.out,c(genus.to.dump,rownames(d.mat)[!mon.gen],add.genus.spec))
	}

#spec.assign	
	if (class(species)!="list") {spec.assign <- spec[assign.id==1]}else{
		spec.assign <- rep(NA,length(spec))
		for (i in 1:length(spec)){
			pom <- species[[assign.id[i]]][i]
			if (is.null(pom)) pom <- NA
			spec.assign[i] <- pom
		}
	}
	pom.tree <- tree.out$tip[!(tree.out$tip %in% spec.assign)]
	for (i in 1:length(spec.assign)){
		if (is.na(spec.assign[i])) { pom <- pom.tree[sel.gen(pom.tree)==sel.gen(spec[i])]
			if (length(pom)>0) spec.assign[i]<- pom
		}
	}
	spec.assign <- cbind(spec,spec.assign)
#randomSpec
	random.trees <- NULL
	if (randomSpec & sum(is.na(spec.assign[,2]))>0){
		random.trees <- list()		
		for (n in 1:ntrees){
			r.spec.assign <- spec.assign
			r.add.genus.spec <- add.genus.spec
			for (i in 1:nrow(r.spec.assign)){
				if (is.na(r.spec.assign[i,2])){
					pom <- r.add.genus.spec[sel.gen(r.add.genus.spec)==sel.gen(r.spec.assign[i,1])]
					if (length(pom)>0) r.spec.assign[i,2]<- sample(pom,1)
					if (!is.na(r.spec.assign[i,2])) r.add.genus.spec <- r.add.genus.spec[r.add.genus.spec!=r.spec.assign[i,2]]
				}
			}
			r.tree.out <- drop.tip(tree.out.add,r.add.genus.spec)
			random.trees <- c(random.trees,list(list(r.tree.out,r.spec.assign)))
		}
	}

	res <- list(subtree=tree.out,spec.assign=spec.assign,random.trees=random.trees)
	return(res)
}

#__________________

#bb <- create_subtree(list(spec.data$Species,spec.data$Acc_species),zanne.tree,randomSpec=T)
#plot(bb[[1]])


