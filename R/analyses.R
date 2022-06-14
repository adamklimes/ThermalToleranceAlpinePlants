# analyses
library(rstan)
library(ape)

# data_loading
dat <- read.csv("data/Dataset_S1.csv")
dat_sp <- read.csv("data/Dataset_S2.csv")
dat_taxonstand <- read.csv("data/analyses/taxonstandNames.csv")
phy <- read.tree("data/phy.tre")

# data_preparation
dat$Species <- with(dat_taxonstand, paste(New.Genus, New.Species)[match(dat$Species, Taxon)])
dat_sp$Species <- with(dat_taxonstand, paste(New.Genus, New.Species)[match(dat_sp$Species, Taxon)])
  # exclusion of duplicated species from another year
dat <- dat[!(dat$Species == "Gentiana algida" & dat$Year == 2014), ]
  # species from 2 localities - preferentially taken from Brennkogel and Antisana
aux <- rowSums(table(dat$Species, dat$Site) > 0)
spec_dup <- names(aux[aux > 1])
dat <- dat[(!dat$Species %in% spec_dup) | dat$Site %in% c("Brennkogel", "Antisana"), ]
dat$Gform <- dat_sp$Growth_form[match(dat$Species, dat_sp$Species)]
dato <- dat[dat$Species %in% names(table(dat$Species)[table(dat$Species) == 4]), ]
# alternative analyses to check robustness of the results
# 1) TTB for all species (including those without 4 obs.)
#   dato <- dat
# 2) TTB with different selection of localities for species from 2 loc.
#   dat <- dat[(!dat$Species %in% spec_dup) | dat$Site %in% c("Yungeno", "Guamaní", "Niwot"), ]
# 3) ad 2 + exclusion of species from "Brennkogel"
#   dat <- dat[dat$Site != "Brennkogel", ]

# data preparation_stan
omax <- tapply(dato$He_Lt_50, dato$Species, max, na.rm = TRUE)
omin <- tapply(dato$Fr_Lt_50, dato$Species, min, na.rm = TRUE)
prep_resp <- function(resp_type, i = 1, dat){
  sp_list <- sort(unique(dat$Species))
  get_period <- function(svar, period){
    sel <- dat$Season == period
    out <- dat[, svar][sel][match(sp_list, dat$Species[sel])]
    names(out) <- sp_list
    out
  }
  if (resp_type == "Fr") return(get_period("Fr_Lt_50", i + 1) - get_period("Fr_Lt_50", i))
  if (resp_type == "He") return(get_period("He_Lt_50", i + 1) - get_period("He_Lt_50", i))
  if (resp_type == "Ra") return(get_period("He_Lt_50", i) - get_period("Fr_Lt_50", i))
}
prep_dat <- function(resp, dat, phy){
  phy$tip.label <- gsub("_", " ", phy$tip.label)
  resp <- resp[!is.na(resp) & is.finite(resp)]
  spec_sel <- names(resp)
  phy_sel <- keep.tip(phy, spec_sel)
  phy_aux <- vcv.phylo(phy_sel)
  phy_vcv <- phy_aux[match(spec_sel, rownames(phy_aux)), match(spec_sel, colnames(phy_aux))]
  dat$Trop <- as.numeric(dat$Site %in% c("Antisana", "Guamaní", "Yungeno"))
  Trop <- tapply(dat$Trop, dat$Species, unique)[spec_sel]
  GrowthF <- as.numeric(factor(tapply(dat$Gform, dat$Species, unique)[spec_sel]))
  Site <- as.numeric(factor(tapply(dat$Site, dat$Species, unique)[spec_sel]))
  Nspec <- length(resp)
  list(Nspec = Nspec, resp = resp, Trop = Trop, GrowthF = GrowthF, Site = Site, phy = phy_vcv)
}
  # dat_stan for each response
dat_stan <- prep_dat(omax-omin, dato, phy)
dat_stan <- prep_dat(omax, dato, phy)
dat_stan <- prep_dat(omin, dato, phy)
dat_stan <- prep_dat(prep_resp("Fr", 1, dat), dat, phy)
dat_stan <- prep_dat(prep_resp("Fr", 2, dat), dat, phy)
dat_stan <- prep_dat(prep_resp("Fr", 3, dat), dat, phy)
dat_stan <- prep_dat(prep_resp("He", 1, dat), dat, phy)
dat_stan <- prep_dat(prep_resp("He", 2, dat), dat, phy)
dat_stan <- prep_dat(prep_resp("He", 3, dat), dat, phy)
  # 4 ranges
dat_stan <- c(prep_dat(prep_resp("Ra", 1, dat), dat, phy),
  prep_dat(prep_resp("Ra", 2, dat), dat, phy),
  prep_dat(prep_resp("Ra", 3, dat), dat, phy),
  prep_dat(prep_resp("Ra", 4, dat), dat, phy))
names(dat_stan) <- paste0(names(dat_stan), rep(1:4, each = 6))

# model fit
  # year_round ranges/max/min
mod_fit <- stan(model_code = mod_code, data = dat_stan, control = list(adapt_delta = 0.99), iter = 4000, seed = 105)
  # differences in consequitive periods
mod_FrHe_fit <- stan(model_code = mod_FrHe_code, data = dat_stan, control = list(adapt_delta = 0.99), iter = 4000, seed = 105)
  # 4 ranges
mod_4ranges_fit <- stan(model_code = mod_4ranges_code, data = dat_stan, iter = 4000, seed = 105)

print(mod_fit, pars ="phylo", include = F)
# save(mod_fit, file = "data/analyses/orange.RData")
# save(mod_fit, file = "data/analyses/omax.RData")
# save(mod_fit, file = "data/analyses/omin.RData")
# save(mod_FrHe_fit, file = "data/analyses/Fr1.RData")
# save(mod_FrHe_fit, file = "data/analyses/Fr2.RData")
# save(mod_FrHe_fit, file = "data/analyses/Fr3.RData")
# save(mod_FrHe_fit, file = "data/analyses/He1.RData")
# save(mod_FrHe_fit, file = "data/analyses/He2.RData")
# save(mod_FrHe_fit, file = "data/analyses/He3.RData")
# save(mod_4ranges_fit, file = "data/analyses/ranges4.RData")

# save(mod_fit, file = "data/analyses/orange_var1.RData")
# save(mod_fit, file = "data/analyses/orange_var2.RData")
# save(mod_fit, file = "data/analyses/orange_var3.RData")

# stan models
mod_code <- "
data {
  int Nspec;
  vector[Nspec] resp;
  vector[Nspec] Trop;
  int GrowthF[Nspec];
  int Site[Nspec];
  cov_matrix[Nspec] phy;
}
parameters {
  real A;
  real B1;
  vector[max(GrowthF)-1] B2;
  real<lower=0> Cpar;
  real<lower=0, upper=1> lambda;
  real<lower=0> sigSite;
  vector[max(Site)] randSite;
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
  sigSite ~ cauchy(0,5);

  for (i in 1:Nspec){
    if (GrowthF[i]>1){pomB2[i] = B2[GrowthF[i]-1];}else{pomB2[i]=0;} 
  }
  randSite ~ normal(0, sigSite);
  pommu = A + B1 * Trop + pomB2 + randSite[Site];
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
  int Site1[Nspec1];
  cov_matrix[Nspec1] phy1;
  int Nspec2;
  vector[Nspec2] resp2;
  vector[Nspec2] Trop2;
  int GrowthF2[Nspec2];
  int Site2[Nspec2];
  cov_matrix[Nspec2] phy2;
  int Nspec3;
  vector[Nspec3] resp3;
  vector[Nspec3] Trop3;
  int GrowthF3[Nspec3];
  int Site3[Nspec3];
  cov_matrix[Nspec3] phy3;
  int Nspec4;
  vector[Nspec4] resp4;
  vector[Nspec4] Trop4;
  int GrowthF4[Nspec4];
  int Site4[Nspec4];
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
  real<lower=0> sigSite;
  vector[max(Site1)] randSite;
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
  sigSite ~ cauchy(0,5);

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
  randSite ~ normal(0, sigSite);
  pommu1 = A1 + B1 * Trop1 + pomC1 + randSite[Site1];
  pommu2 = A2 + B2 * Trop2 + pomC2 + randSite[Site2];
  pommu3 = A3 + B3 * Trop3 + pomC3 + randSite[Site3];
  pommu4 = A4 + B4 * Trop4 + pomC4 + randSite[Site4];
  resp1 ~ multi_normal(pommu1, phylo1);
  resp2 ~ multi_normal(pommu2, phylo2);
  resp3 ~ multi_normal(pommu3, phylo3);
  resp4 ~ multi_normal(pommu4, phylo4);
}
"
#_
