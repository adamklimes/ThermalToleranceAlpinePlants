# Analysis_latitude_range
library(caper)

# loading data
rapdat <- read.csv("data/Rapoport_dat.csv")
load(file = "data/analyses/r_tree.RData")

# data preparation
rapdat$Spec <- gsub(" ", "_", rapdat$Spec)
rapdat <- rapdat[!is.na(rapdat$Range_resist), ]
tree[[1]]$node.label <- NULL

# species corrections
rapdat$Spec[rapdat$Spec == "Cerastium_arvense_s.l."] <- "Cerastium_arvense"
rapdat$Spec[rapdat$Spec == "Cerastium_uniflorum"] <- "Cerastium_cf."
rapdat$Spec[rapdat$Spec == "Vaccinium_cf._caespitosum"] <- "Vaccinium_cf."
rapdat$Spec[rapdat$Spec == "Senecio_hohenackeri_=_\"lasioovatus\""] <- "Senecio_lasioovatus2"
rapdat$Spec[rapdat$Spec == "Senecio_rhizomatus_=_\"formosus\""] <- "Senecio_formosus2"
rapdat$Spec[rapdat$Spec == "Valeriana_micropterina_=_\"compuesta\""] <- "Valeriana_compuesta2"
rapdat$Spec[rapdat$Spec == "Loricaria_antisanensis"] <- "Gamochaeta_sp"
rapdat$Spec[rapdat$Spec == "Oritrophium_hieracioides"] <- "Diplostephium_sp2"
rapdat$Spec[rapdat$Spec == "Oritrophium_peruvianum"] <- "Diplostephium_sp"
rapdat$Spec[rapdat$Spec == "Lasiocephalus_ovatus"] <- "Senecio_sp"
rapdat$Spec[rapdat$Spec == "Pentacalia_andicola"] <- "Senecio_sp2"
rapdat$Spec[rapdat$Spec == "Silene_\"rimbachii\""] <- "Silene_cf."
rapdat$Spec[rapdat$Spec == "Silene_cf._chilensis"] <- "Silene_thysanodes2"

# assigned species names
rapdat$Spec_assign <- tree[[2]][match(rapdat$Spec, tree[[2]][, 1]), 2]

# analysis
dat <- comparative.data(
  phy = tree[[1]],
  data = rapdat,
  names.col = Spec_assign,
  vcv = TRUE, 
  na.omit = FALSE 
  # ,warn.dropped = TRUE # not all species are in rapdat 
)

mod_e <- pgls(sqrt(Range_equator) ~ Tropical * Range_resist, data = dat, 
  lambda = "ML")

mod_l <- pgls(sqrt(Range_latitude) ~ Tropical * Range_resist, data = dat, 
  lambda = "ML")

# sensitivity to one observation with large resist
dat2 <- comparative.data(
  phy = tree[[1]],
  data = rapdat[rapdat$Range_resist < 90, ],
  names.col = Spec_assign,
  vcv = TRUE, 
  na.omit = FALSE 
  # ,warn.dropped = TRUE # not all species are in rapdat 
)

mod_l2 <- pgls(sqrt(Range_latitude) ~ Tropical * Range_resist, data = dat2, 
  lambda = "ML")

# png("figures/Fig5_resist_lat.png", width = 480*10, height = 480*10, res = 72*10)
cols <- -rapdat$Tropical*2+4
par(mai = c(0.9,0.8,0.1,0.1))
plot(sqrt(Range_latitude) ~ Range_resist, data = rapdat, col = cols, 
  pch = c(16,1)[(Range_resist > 90) + 1], axes = FALSE, 
  xlab = "Temperature resistance range [°C]", ylab = "Latitudinal range [°]")
with(rapdat[rapdat$Range_resist > 90, ],
  text(Range_resist, sqrt(Range_latitude), gsub("_", " ", Spec), pos = 2, 
  cex = 0.8))
xx_s <- with(rapdat[rapdat$Tropical == 0, ], 
  seq(min(Range_resist), max(Range_resist), length.out = 100))
xx_t <- with(rapdat[rapdat$Tropical == 1, ], 
  seq(min(Range_resist), max(Range_resist), length.out = 100))
lines(xx_s, coef(mod_l)[1] + coef(mod_l)[3] * xx_s, col = cols[1], lwd = 2)
lines(xx_t, sum(coef(mod_l)[c(1,2)]) + sum(coef(mod_l)[c(3,4)]) * xx_t, col = cols[80],
  lwd = 2)
xx_s2 <- with(rapdat[rapdat$Tropical == 0 & rapdat$Range_resist < 90, ], 
  seq(min(Range_resist), max(Range_resist), length.out = 100))
lines(xx_s2, coef(mod_l2)[1] + coef(mod_l2)[3] * xx_s2, lty = 2, lwd = 2, col = cols[1])
box(bty = "l")
axis(1)
axis(2, labels = c(1,5,10,20,30,50), at = sqrt(c(1,5,10,20,30,50)))
legend("bottomright", pch = 16, col = cols[c(1,80)], 
  legend = c("Temperate", "Tropical"), bty = "n")

#_