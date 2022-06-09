#stan_model
library(rstan)
library(ape)

# loading data
dat <- read.table("data/Temperature_resistance_data.csv", sep = ";", header = TRUE)
spec_gform <- read.table("data/Growth_Forms_List.csv", sep = ";", header = TRUE)
load(file = "data/analyses/r_tree.RData")

# data preparation
dat$gform <- as.character(spec_gform$Growth_form[match(dat$Species, spec_gform$Species)])
  # exclusion of duplicated species from another year
dat <- dat[!(dat$Species == "Gentiana_algida" & dat$Year == 2014), ]
  # species names
dat$Species[dat$Species == "Pedicularis_parryi_late_melt"] <- "Pedicularis_parryi2"
pom_names <- unlist(lapply(strsplit(dat$Species, "_"), function(x) paste(x[1], x[2], sep = "_")))
pom_names[grepl("\"", pom_names)] <- paste0(pom_names[grepl("\"", pom_names)], 2)
pom_names <- gsub("\"", "", pom_names)
old_spec <- c("Loricaria_antisanensis", "Oritrophium_peruvianum", "Oritrophium_hieracioides", "Lasiocephalus_ovatus", "Pentacalia_andicola")
new_spec <- c("Gamochaeta_sp", "Diplostephium_sp", "Diplostephium_sp2", "Senecio_sp", "Senecio_sp2")
pom_names[pom_names %in% old_spec] <- new_spec[match(pom_names[pom_names %in% old_spec], old_spec)]
dat$Species <- pom_names

  # selected species assesed in 4 periods
  # species from 2 localities - randomly chosen one of the localities 
sel_spec_pom <- colSums(table(dat$Period, dat$Species) > 0)
sel_spec <- names(sel_spec_pom)[sel_spec_pom > 3]
dat_cons <- dat[(dat$Species %in% sel_spec) & !(is.na(dat$Fr_lt50) & is.na(dat$He_lt50)), ]
set.seed(10)
for (i in unique(dat_cons$Species)){
  loc <- sample(unique(dat_cons$Site[dat_cons$Species == i]), 1)
  dat_cons <- dat_cons[dat_cons$Species != i | dat_cons$Site == loc, ]  
}

# data preparation_stan
omax <- tapply(dat_cons$He_lt50, dat_cons$Species, max, na.rm = TRUE)
omin <- tapply(dat_cons$Fr_lt50, dat_cons$Species, min, na.rm = TRUE)
prep_resp <- function(resp_type, i = 1, dat_cons){
  sp_list <- sort(unique(dat_cons$Species))
  get_period <- function(svar, period){
    sel <- dat_cons$Period == period
    out <- dat_cons[, svar][sel][match(sp_list, dat_cons$Species[sel])]
    names(out) <- sp_list
    out
  }
  if (resp_type == "Fr") return(get_period("Fr_lt50", i + 1) - get_period("Fr_lt50", i))
  if (resp_type == "He") return(get_period("He_lt50", i + 1) - get_period("He_lt50", i))
  if (resp_type == "Ra") return(get_period("He_lt50", i) - get_period("Fr_lt50", i))
}
prep_dat <- function(resp, dat_cons, tree){
  spec_aval <- !(resp == -Inf | resp == Inf | is.na(resp))
  resp <- resp[spec_aval]
  pres_spec <- tree[[2]][tree[[2]][, 1] %in% dat_cons$Species, ]
  cons_r_tree <- drop.tip(tree[[1]], tree[[1]]$tip.label[!(tree[[1]]$tip.label %in% pres_spec[, 2])])
  phy_pom <- vcv.phylo(cons_r_tree)
  phy <- phy_pom[match(pres_spec[, 2], colnames(phy_pom))[order(pres_spec[, 1])], match(pres_spec[, 2], colnames(phy_pom))[order(pres_spec[, 1])]]
  phy <- phy[spec_aval, spec_aval]
  dat_cons$Trop <- as.numeric(dat_cons$Site %in% c("ANT", "GUA", "Yungeno"))
  Trop <- tapply(dat_cons$Trop, dat_cons$Species, mean)[spec_aval]
  GrowthF <- as.numeric(factor(tapply(dat_cons$gform, dat_cons$Species, unique)[spec_aval]))
  Area <- as.numeric(factor(tapply(dat_cons$Site, dat_cons$Species, unique)[spec_aval]))
  Nspec <- length(Area)
  list(Nspec = Nspec, resp = resp, Trop = Trop, GrowthF = GrowthF, Area = Area, phy = phy)
}
  # dat_stan for each response
dat_stan <- prep_dat(omax-omin, dat_cons, tree)
dat_stan <- prep_dat(omax, dat_cons, tree)
dat_stan <- prep_dat(omin, dat_cons, tree)
dat_stan <- prep_dat(prep_resp("Fr", 1, dat_cons), dat_cons, tree)
dat_stan <- prep_dat(prep_resp("Fr", 2, dat_cons), dat_cons, tree)
dat_stan <- prep_dat(prep_resp("Fr", 3, dat_cons), dat_cons, tree)
dat_stan <- prep_dat(prep_resp("He", 1, dat_cons), dat_cons, tree)
dat_stan <- prep_dat(prep_resp("He", 2, dat_cons), dat_cons, tree)
dat_stan <- prep_dat(prep_resp("He", 3, dat_cons), dat_cons, tree)
  # 4 ranges
dat_stan <- c(prep_dat(prep_resp("Ra", 1, dat_cons), dat_cons, tree),
  prep_dat(prep_resp("Ra", 2, dat_cons), dat_cons, tree),
  prep_dat(prep_resp("Ra", 3, dat_cons), dat_cons, tree),
  prep_dat(prep_resp("Ra", 4, dat_cons), dat_cons, tree))
names(dat_stan) <- paste0(names(dat_stan), rep(1:4, each = 6))

# model fit
  # year_round ranges/max/min
mod_fit <- stan(model_code = mod_code, data = dat_stan, control = list(adapt_delta = 0.99), iter = 4000, seed = 105)
  # differences in consequitive periods
mod_FrHe_fit <- stan(model_code = mod_FrHe_code, data = dat_stan, control = list(adapt_delta = 0.99), iter = 4000, seed = 105)
  # 4 ranges
mod_4ranges_fit <- stan(model_code = mod_4ranges_code, data = dat_stan_ranges, iter = 4000, seed = 105)

# stan models
mod_code <- "
data {
  int Nspec;
  vector[Nspec] resp;
  vector[Nspec] Trop;
  int GrowthF[Nspec];
  int Area[Nspec];
  cov_matrix[Nspec] phy;
}
parameters {
  real A;
  real B1;
  vector[max(GrowthF)-1] B2;
  real<lower=0> Cpar;
  real<lower=0, upper=1> lambda;
  real<lower=0> sigArea;
  vector[max(Area)] randArea;
}
transformed parameters{
  matrix[Nspec, Nspec] phylo;
  phylo = Cpar*((phy-diag_matrix(diagonal(phy)))*lambda + diag_matrix(diagonal(phy)));
}
model {
  vector[Nspec] pommu;
  vector[Nspec] pomB2;
  //priors
  B1 ~ cauchy(0,5);
  B2 ~ cauchy(0,5);
  Cpar ~ cauchy(0,5);
  sigArea ~ cauchy(0,5);

  for (i in 1:Nspec){
    if (GrowthF[i]>1){pomB2[i] = B2[GrowthF[i]-1];}else{pomB2[i]=0;} 
  }
  randArea ~ normal(0, sigArea);
  pommu = A + B1 * Trop + pomB2 + randArea[Area];
  resp ~ multi_normal(pommu, phylo);
}
"

mod_FrHe_code <- "
data {
  int Nspec;
  vector[Nspec] resp;
  vector[Nspec] Trop;
  cov_matrix[Nspec] phy;
}
parameters {
  real A;
  real B;
  real<lower=0> Cpar;
  real<lower=0, upper=1> lambda;
}
transformed parameters{
  matrix[Nspec, Nspec] phylo;
  phylo = Cpar*((phy-diag_matrix(diagonal(phy)))*lambda + diag_matrix(diagonal(phy)));
}
model {
  vector[Nspec] pommu;
  //priors
  B ~ cauchy(0,5);
  Cpar ~ cauchy(0,5);

  pommu = A + B * Trop; 
  resp ~ multi_normal(pommu, phylo);
}
"

mod_4ranges_code <- "
data {
  int Nspec1;
  vector[Nspec1] resp1;
  vector[Nspec1] Trop1;
  int GrowthF1[Nspec1];
  int Area1[Nspec1];
  cov_matrix[Nspec1] phy1;
  int Nspec2;
  vector[Nspec2] resp2;
  vector[Nspec2] Trop2;
  int GrowthF2[Nspec2];
  int Area2[Nspec2];
  cov_matrix[Nspec2] phy2;
  int Nspec3;
  vector[Nspec3] resp3;
  vector[Nspec3] Trop3;
  int GrowthF3[Nspec3];
  int Area3[Nspec3];
  cov_matrix[Nspec3] phy3;
  int Nspec4;
  vector[Nspec4] resp4;
  vector[Nspec4] Trop4;
  int GrowthF4[Nspec4];
  int Area4[Nspec4];
  cov_matrix[Nspec4] phy4;
}
parameters {
  real A1;
  real A2;
  real A3;
  real A4;
  real B1;
  real B2;
  real B3;
  real B4;
  vector[max(GrowthF1)-1] C;
  real<lower=0> Cpar1;
  real<lower=0> Cpar2;
  real<lower=0> Cpar3;
  real<lower=0> Cpar4;
  real<lower=0, upper=1> lambda1;
  real<lower=0, upper=1> lambda2;
  real<lower=0, upper=1> lambda3;
  real<lower=0, upper=1> lambda4;
  real<lower=0> sigArea;
  vector[max(Area1)] randArea;
}
model {
  vector[Nspec1] pommu1;
  vector[Nspec2] pommu2;
  vector[Nspec3] pommu3;
  vector[Nspec4] pommu4;
  vector[Nspec1] pomC1;
  vector[Nspec2] pomC2;
  vector[Nspec3] pomC3;
  vector[Nspec4] pomC4;
  matrix[Nspec1, Nspec1] phylo1;
  matrix[Nspec2, Nspec2] phylo2;
  matrix[Nspec3, Nspec3] phylo3;
  matrix[Nspec4, Nspec4] phylo4;
  //priors
  B1 ~ cauchy(0,5);
  B2 ~ cauchy(0,5);
  B3 ~ cauchy(0,5);
  B4 ~ cauchy(0,5);
  C ~ cauchy(0,5);
  Cpar1 ~ cauchy(0,5);
  Cpar2 ~ cauchy(0,5);
  Cpar3 ~ cauchy(0,5);
  Cpar4 ~ cauchy(0,5);
  sigArea ~ cauchy(0,5);

  phylo1 = Cpar1*((phy1-diag_matrix(diagonal(phy1)))*lambda1 + diag_matrix(diagonal(phy1)));
  phylo2 = Cpar2*((phy2-diag_matrix(diagonal(phy2)))*lambda2 + diag_matrix(diagonal(phy2)));
  phylo3 = Cpar3*((phy3-diag_matrix(diagonal(phy3)))*lambda3 + diag_matrix(diagonal(phy3)));
  phylo4 = Cpar4*((phy4-diag_matrix(diagonal(phy4)))*lambda4 + diag_matrix(diagonal(phy4)));

  for (i in 1:Nspec1){
    if (GrowthF1[i]>1){pomC1[i] = C[GrowthF1[i]-1];}else{pomC1[i]=0;} 
  }
  for (i in 1:Nspec2){
    if (GrowthF2[i]>1){pomC2[i] = C[GrowthF2[i]-1];}else{pomC2[i]=0;} 
  }
  for (i in 1:Nspec3){
    if (GrowthF3[i]>1){pomC3[i] = C[GrowthF3[i]-1];}else{pomC3[i]=0;} 
  }
  for (i in 1:Nspec4){
    if (GrowthF4[i]>1){pomC4[i] = C[GrowthF4[i]-1];}else{pomC4[i]=0;} 
  }
  randArea ~ normal(0, sigArea);
  pommu1 = A1 + B1 * Trop1 + pomC1 + randArea[Area1];
  pommu2 = A2 + B2 * Trop2 + pomC2 + randArea[Area2];
  pommu3 = A3 + B3 * Trop3 + pomC3 + randArea[Area3];
  pommu4 = A4 + B4 * Trop4 + pomC4 + randArea[Area4];
  resp1 ~ multi_normal(pommu1, phylo1);
  resp2 ~ multi_normal(pommu2, phylo2);
  resp3 ~ multi_normal(pommu3, phylo3);
  resp4 ~ multi_normal(pommu4, phylo4);
}
"

#_
