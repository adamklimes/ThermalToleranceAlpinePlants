# Analysis_latitude_range
library(caper)

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
phy$node.label <- NULL
omax <- tapply(dato$He_Lt_50, dato$Species, max, na.rm = TRUE)
omin <- tapply(dato$Fr_Lt_50, dato$Species, min, na.rm = TRUE)
orange <- omax - omin
numcoor <- function(x) {
  hemis <- gsub("[[:digit:]]+\\.*", "", x)
  cf <- ifelse(hemis == "S", -1, 1)
  cf * as.numeric(sub("[[:upper:]]", "", x))
}
aux <- numcoor(dat_sp$Northern_limit)
dat_sp$Range_latitude <- aux - numcoor(dat_sp$Southern_limit)
dat_sp$Range_equator <- aux - ifelse(aux > 0, pmax(0, numcoor(dat_sp$Southern_limit)), numcoor(dat_sp$Southern_limit))
dat_sp$orange <- orange[match(dat_sp$Species, names(orange))]
lmax <- tapply(dat$He_Lt_50, dat$Species, max, na.rm = TRUE)
lmin <- tapply(dat$Fr_Lt_50, dat$Species, min, na.rm = TRUE)
dat_sp$lrange <- (lmax - lmin)[match(dat_sp$Species, names(lmax))]
dat_sp$lrange[!is.finite(dat_sp$lrange)] <- NA
dat_sp$Species <- gsub(" ", "_", dat_sp$Species)
dat_sp$Tropical <- as.numeric(dat_sp$Locality %in% c("EC", "BOL", "EC_BOL"))
# tapply(dat_sp$Range_latitude, dat_sp$Tropical, mean, na.rm = TRUE)
# analysis
cdat <- comparative.data(
  phy = phy,
  data = dat_sp,
  names.col = Species,
  vcv = TRUE, 
  na.omit = FALSE 
)

mod_e <- pgls(sqrt(Range_equator) ~ Tropical * orange, data = cdat, 
  lambda = "ML")
mod_l <- pgls(sqrt(Range_latitude) ~ Tropical * orange, data = cdat, 
  lambda = "ML")
mod_le <- pgls(sqrt(Range_equator) ~ Tropical * lrange, data = cdat, 
  lambda = "ML")
mod_ll <- pgls(sqrt(Range_latitude) ~ Tropical * lrange, data = cdat, 
  lambda = "ML")
mod_l2 <- pgls(sqrt(Range_latitude) ~ Tropical * orange, data = cdat[-71,], 
  lambda = "ML") # without Cerastium uniflorum
anova(mod_e)
anova(mod_l)
anova(mod_le)
anova(mod_ll)
anova(mod_l2)

aux_line <- function(x, mod, trop, col, lty = 1){
  xx <- seq(min(dat_sp[dat_sp$Tropical == trop, x], na.rm = TRUE), 
    max(dat_sp[dat_sp$Tropical == trop, x], na.rm = TRUE), length.out = 100)
  if (trop) lines(xx, sum(coef(mod)[c(1,2)]) + sum(coef(mod)[c(3,4)]) * xx, 
    col = col, lwd = 2, lty = lty) else
    lines(xx, coef(mod)[1] + coef(mod)[3] * xx, col = col, lwd = 2, lty = lty)
}

# png("figures/Fig5_resist_lat.png", width = 480*10, height = 480*10, res = 72*10)
cols <- c("blue", "red")
par(mai = c(0.9,0.8,0.1,0.1))
plot(sqrt(Range_latitude) ~ lrange, data = dat_sp, col = cols[dat_sp$Tropical + 1], 
  pch = c(16, 1)[is.na(dat_sp$orange) + 1], axes = FALSE, 
  xlab = "Thermal tolerance breadth [°C]", ylab = "Latitudinal range [°]")
aux_line("orange", mod_l, trop = 0, col = cols[1])
aux_line("orange", mod_l, trop = 1, col = cols[2])
aux_line("lrange", mod_ll, trop = 0, col = cols[1], lty = 2)
aux_line("lrange", mod_ll, trop = 1, col = cols[2], lty = 2)
box(bty = "l")
axis(1)
axis(2, labels = c(1,5,10,20,30,50), at = sqrt(c(1,5,10,20,30,50)))
legend("bottomright", pch = 16, col = cols, lwd = 2,
  legend = c("Temperate", "Tropical"), bty = "n")

#_