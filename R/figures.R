# figures
# follows script "analyses.R"

#data preparation
dat_cons$Trop <- as.numeric(dat_cons$Site %in% c("ANT", "GUA", "Yungeno"))
dat_cons$Site <- factor(dat_cons$Site, levels = unique(dat_cons$Site)[c(1,5,2:4)])
dat_cons$gform <- factor(dat_cons$gform, levels = unique(dat_cons$gform)[c(1:2,4,3,5:6)])
cols <- c("blue", "red")

## Figures
# png("figures/Fig3_FrHe.png", height = 480*10, width = 480*10, res = 72*10)
par(mai = c(1,0.8,0.1,0.1))
boxplot(Fr_lt50 ~ Trop * Period, dat_cons, col = cols,
  xlab = "Period", ylab = "Temperature [°C]", ylim = c(-21, 25), 
  at = c(2,3,5,6,8,9,11,12)/2, axes = FALSE)
boxplot((He_lt50-30) ~ Trop * Period, dat_cons, col = cols, 
  add = TRUE, at = c(2,3,5,6,8,9,11,12)/2, axes = FALSE)
axis(2, labels = c(-50,-20,-15,-10,-5,0), at = c(-50,-20,-15,-10,-5,0))
axis(2, labels = c(35,40,45,50,55), at = c(5,10,15,20,25))
axis(1, at = 1:4*1.5-0.25, labels = 1:4, tick = FALSE)
axis(1, at = c(-100,100), labels = c("", ""))
sq <- seq(0,5, length.out = 8)
axis(2, at = sq[2:3], labels = c("",""), lwd.ticks = -1)
axis(2, at = sq[4:5], labels = c("",""), lwd.ticks = -1)
axis(2, at = sq[6:7], labels = c("",""), lwd.ticks = -1)
abline(h = mean(sq[4:5]), lty = 2)

#____________________________________________________
# Fig 2 - Results
library(rstan)
load(file = "data/analyses/Fr1.RData")
post_Fr1 <- extract(mod_FrHe_fit, pars = c("phylo","Cpar","lp__"), include = FALSE)
load(file = "data/analyses/Fr2.RData")
post_Fr2 <- extract(mod_FrHe_fit, pars = c("phylo","Cpar","lp__"), include = FALSE)
load(file = "data/analyses/Fr3.RData")
post_Fr3 <- extract(mod_FrHe_fit, pars = c("phylo","Cpar","lp__"), include = FALSE)
load(file = "data/analyses/He1.RData")
post_He1 <- extract(mod_FrHe_fit, pars = c("phylo","Cpar","lp__"), include = FALSE)
load(file = "data/analyses/He2.RData")
post_He2 <- extract(mod_FrHe_fit, pars = c("phylo","Cpar","lp__"), include = FALSE)
load(file = "data/analyses/He3.RData")
post_He3 <- extract(mod_FrHe_fit, pars = c("phylo","Cpar","lp__"), include = FALSE)
load(file = "data/analyses/orange.RData")
post_range <- extract(mod_fit, pars = c("phylo","Cpar","lp__"), include = FALSE)
load(file = "data/analyses/ranges4.RData")
post_ranges4 <- extract(mod_4ranges_fit, pars = c("phylo","Cpar","lp__"), include = FALSE)

aux_lines <- function(x, mod, y_shift = 0, lwds = 1:3){
  x[2] <- x[2] - 0.1
  lines(x-c(0,0.05), c(0,mean(mod$A)) + y_shift, lwd = lwds[2], col = cols[1])
  lines(x+c(0,0.05), c(0,mean(mod$A + mod$B)) + y_shift, col = cols[2], lwd = lwds[2])
  arrows(rep(x[2], 2)+c(-0.05,0.05), c(quantile(mod$A, 0.025), 
    quantile(mod$A + mod$B, 0.025)) + y_shift, y1 = c(quantile(mod$A, 0.975), 
    quantile(mod$A + mod$B, 0.975)) + y_shift, code = 3, length = 0, col = cols,
    lwd = lwds[1])
  arrows(rep(x[2], 2) + c(-0.05,0.05), c(quantile(mod$A, 0.25), 
    quantile(mod$A + mod$B, 0.25)) + y_shift, y1 = c(quantile(mod$A, 0.75), 
    quantile(mod$A + mod$B, 0.75)) + y_shift, code = 3, length = 0, col = cols, 
    lwd = lwds[3])
}
temp_range <- with(post_range, A + rowSums(B2) / 6)
trop_range <- with(post_range, A + B1 + rowSums(B2) / 6)
temp_range1 <- with(post_ranges4, A1 + rowSums(C) / 6)
trop_range1 <- with(post_ranges4, A1 + B1 + rowSums(C) / 6)
temp_range2 <- with(post_ranges4, A2 + rowSums(C) / 6)
trop_range2 <- with(post_ranges4, A2 + B2 + rowSums(C) / 6)
temp_range3 <- with(post_ranges4, A3 + rowSums(C) / 6)
trop_range3 <- with(post_ranges4, A3 + B3 + rowSums(C) / 6)
temp_range4 <- with(post_ranges4, A4 + rowSums(C) / 6)
trop_range4 <- with(post_ranges4, A4 + B4 + rowSums(C) / 6)
temp_ranges <- cbind(temp_range1, temp_range2, temp_range3, temp_range4)
trop_ranges <- cbind(trop_range1, trop_range2, trop_range3, trop_range4)

# png("figures/Fig4_Results.png", height = 480*10, width = 480*10, res = 72*10)
layout(matrix(c(1,1,2,3,2,3), 2, 3))
par(mai = c(0.1,0.8,0.1,0.1))
plot(1:2, c(0,80), type = "n", axes = FALSE, 
  ylab = "Thermal tolerance breadth [K]", xlab = "", cex.lab = 1.3)
axis(2, labels = 0:8*10, at = 0:8*10, cex.axis = 1.3)
axis(2, labels = c("",""), at = c(0,100))
barplot(c(mean(temp_range), mean(trop_range)), add = TRUE, col = cols, width = 0.5, 
  space = c(2,0), axes = FALSE)
arrows(1.25, quantile(temp_range, 0.025), 
  y1 = quantile(temp_range, 0.975), code = 3, angle = 90, length = 0.12)
arrows(1.25, quantile(temp_range, 0.25), 
  y1 = quantile(temp_range, 0.75), code = 3, length = 0, lwd = 3)
arrows(1.75, quantile(trop_range, 0.025), 
  y1 = quantile(trop_range, 0.975), code = 3, angle = 90, length = 0.12)
arrows(1.75, quantile(trop_range, 0.25), 
  y1 = quantile(trop_range, 0.75), code = 3, length = 0, lwd = 3)
par(mai = c(0.43,0.8,0.1,0.1))
plot(c(0.5,6), c(46,68), type = "n", axes = FALSE, 
  ylab = "Thermal tolerance breadth [K]", xlab = "", cex.lab = 1.3)
axis(2, labels = 9:14*5, at = 9:14*5, cex.axis = 1.3)
axis(2, labels = c("",""), at = c(48,100), tick = FALSE)
barplot(rbind(colMeans(temp_ranges), colMeans(trop_ranges)) - 46, 
  add = TRUE, col = cols, width = 0.5, 
  axes = FALSE, beside = TRUE, offset = 46, names.arg = rep("",8))
arrows(1:4*1.5-0.75, apply(temp_ranges, 2, quantile, 0.025), 
  y1 = apply(temp_ranges, 2, quantile, 0.925), 
  code = 3, angle = 90, length = 0.07)
arrows(1:4*1.5-0.75, apply(temp_ranges, 2, quantile, 0.25), 
  y1 = apply(temp_ranges, 2, quantile, 0.75), 
  code = 3, length = 0, lwd = 3)
arrows(1:4*1.5-0.25, apply(trop_ranges, 2, quantile, 0.025), 
  y1 = apply(trop_ranges, 2, quantile, 0.925), 
  code = 3, angle = 90, length = 0.07)
arrows(1:4*1.5-0.25, apply(trop_ranges, 2, quantile, 0.25), 
  y1 = apply(trop_ranges, 2, quantile, 0.75), 
  code = 3, length = 0, lwd = 3)
text(1:4*1.5-0.5, 67, paste("Period", 1:4), cex = 1.3)
par(mai = c(0.6,0.8,0.3,0.1))
plot(c(1,4), c(-8,18), type="n", xlab="Period", ylab="", axes = FALSE, 
  cex.lab = 1.3)
abline(h = c(0, 12), lty = 2, col = "grey")
aux_lines(1:2, post_Fr1)
aux_lines(2:3, post_Fr2)
aux_lines(3:4, post_Fr3)
aux_lines(1:2, post_He1, y_shift = 12)
aux_lines(2:3, post_He2, y_shift = 12)
aux_lines(3:4, post_He3, y_shift = 12)
axis(1, labels = 1:4, at = 1:4, cex.axis = 1.3)
axis(1, labels = c("", ""), at = c(-10, 10))
sq <- seq(5,8, length.out = 8)
axis(2, at = sq[2:3], labels = c("",""), lwd.ticks = -1)
axis(2, at = sq[4:5], labels = c("",""), lwd.ticks = -1)
axis(2, at = sq[6:7], labels = c("",""), lwd.ticks = -1)
abline(h = mean(sq[4:5]), lty = 2)
axis(2, labels = c(-3,0,3), at = c(-3,0,3), cex.axis = 1.3)
axis(2, labels = c(-3,0,3), at = c(-3,0,3) + 12, cex.axis = 1.3)
axis(2, labels = c("", ""), at = c(-10, sq[1] - 0.2), lwd.ticks = -1)
axis(2, labels = c("", ""), at = c(sq[8] + 0.2, 100), lwd.ticks = -1)
axis(2, labels = c("Heat", "Freezing"), at = c(12, 0), tick = FALSE, 
  line = 1.4, cex.axis = 1.3)
axis(2, labels = expression(paste(Delta, "K")), at = mean(sq[4:5]), 
  tick = FALSE, line = 2.4, cex.axis = 1.3)
par(new = TRUE, mai = c(0,0,0,0), mfrow = c(1,1))
plot(0:1,0:1, axes = FALSE, xaxs = "i", yaxs = "i", ann = FALSE, type = "n")
text(c(0.04,0.36,0.36), c(0.95,0.95,0.5), paste0("H",1:3,")"), cex = 1.6)
legend(0.45,0.57, fill=cols, legend=c("Temperate","Tropics"), bty = "n")

#____________________________________________________
## Figure 4 - components
# H1 + H2
off <- 45
# png("figures/Fig4_H1H2.png", height = 480*10, width = 480*10*2, res = 72*10)
par(mai = c(0.1,0.8,0.1,0.1))
plot(c(-1.75,6), c(off,75), type = "n", axes = FALSE, 
  ylab = "", xlab = "")
axis(2, labels = 0:16*5, at = 0:16*5, cex.axis = 2.5, lwd = 11)
axis(2, labels = c("",""), at = c(0,100), lwd = 11)
barplot(c(mean(temp_range), mean(trop_range)) - off, add = TRUE, col = cols, 
  width = 0.5, space = c(-3,0), axes = FALSE, offset = off)
arrows(-1.25, quantile(temp_range, 0.025), lwd = 5,
  y1 = quantile(temp_range, 0.975), code = 3, angle = 90, length = 0.07)
arrows(-1.25, quantile(temp_range, 0.25), 
  y1 = quantile(temp_range, 0.75), code = 3, length = 0, lwd = 11)
arrows(-0.75, quantile(trop_range, 0.025), lwd = 5,
  y1 = quantile(trop_range, 0.975), code = 3, angle = 90, length = 0.07)
arrows(-0.75, quantile(trop_range, 0.25), 
  y1 = quantile(trop_range, 0.75), code = 3, length = 0, lwd = 11)
barplot(rbind(colMeans(temp_ranges), colMeans(trop_ranges)) - off, 
  add = TRUE, col = cols, width = 0.5, offset = off,
  axes = FALSE, beside = TRUE, names.arg = rep("",8))
arrows(1:4*1.5-0.75, apply(temp_ranges, 2, quantile, 0.025), 
  y1 = apply(temp_ranges, 2, quantile, 0.925), 
  code = 3, angle = 90, length = 0.07, lwd = 5)
arrows(1:4*1.5-0.75, apply(temp_ranges, 2, quantile, 0.25), 
  y1 = apply(temp_ranges, 2, quantile, 0.75), 
  code = 3, length = 0, lwd = 11)
arrows(1:4*1.5-0.25, apply(trop_ranges, 2, quantile, 0.025), 
  y1 = apply(trop_ranges, 2, quantile, 0.925), 
  code = 3, angle = 90, length = 0.07, lwd = 5)
arrows(1:4*1.5-0.25, apply(trop_ranges, 2, quantile, 0.25), 
  y1 = apply(trop_ranges, 2, quantile, 0.75), 
  code = 3, length = 0, lwd = 11)

#________________________
# H3
# png("figures/Fig4_H3.png", height = 480*10, width = 480*10*2, res = 72*10)
par(mai = c(0.1,0.8,0.1,0.1))
plot(c(-0.3,4.1), c(-8,18), type="n", xlab="", ylab="", axes = FALSE)
#abline(h = c(0, 12), lty = 2, col = "grey")
lines(c(0.65,4.15), c(0,0), lwd = 3, lty = 2)
lines(c(0.65,4.15), c(12,12), lwd = 3, lty = 2)
aux_lines(1:2, post_Fr1, lwds = c(5,8,11))
aux_lines(2:3, post_Fr2, lwds = c(5,8,11))
aux_lines(3:4, post_Fr3, lwds = c(5,8,11))
aux_lines(1:2, post_He1, y_shift = 12, lwds = c(5,8,11))
aux_lines(2:3, post_He2, y_shift = 12, lwds = c(5,8,11))
aux_lines(3:4, post_He3, y_shift = 12, lwds = c(5,8,11))
#axis(1, labels = 1:4, at = 1:4, cex.axis = 1.3)
#axis(1, labels = c("", ""), at = c(-10, 10))
sq <- seq(5,8, length.out = 8)
axis(2, at = sq[2:3], labels = c("",""), lwd.ticks = -1, lwd = 11)
axis(2, at = sq[4:5], labels = c("",""), lwd.ticks = -1, lwd = 11)
axis(2, at = sq[6:7], labels = c("",""), lwd.ticks = -1, lwd = 11)
abline(h = mean(sq[4:5]), lty = 2, lwd = 5)
axis(2, labels = c(-3,0,3), at = c(-3,0,3), cex.axis = 2.5, lwd = 11)
axis(2, labels = c(-3,0,3), at = c(-3,0,3) + 12, cex.axis = 2.5, lwd = 11)
axis(2, labels = c("", ""), at = c(-10, sq[1] - 0.2), lwd.ticks = -1, lwd = 11)
axis(2, labels = c("", ""), at = c(sq[8] + 0.2, 100), lwd.ticks = -1, lwd = 11)

#____________________________________________________
# growth forms
get_ci <- function(x){
  q <- quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975))
  q[3] <- mean(x)
  names(q)[3] <- "mean"
  q
}
load(file="data/analyses/omax.RData")
post_max <- extract(mod_fit, pars=c("phylo","Cpar","lp__"), include=F)
load(file="data/analyses/omin.RData")
post_min <- extract(mod_fit, pars=c("phylo","Cpar","lp__"), include=F)

gf_range <- with(post_range, as.vector(A) + cbind(0, B2) + as.vector(B1/2))
gf_max <- with(post_max, as.vector(A) + cbind(0, B2) + as.vector(B1/2))
gf_min <- with(post_min, as.vector(A) + cbind(0, B2) + as.vector(B1/2))
gfra <- apply(gf_range, 2, get_ci)
gfma <- apply(gf_max, 2, get_ci)
gfmi <- apply(gf_min, 2, get_ci)
plot_ci <- function(mat, x){
  points(x, mat["mean", ], pch = 16, cex = 1.3)
  arrows(x, mat[1, ], y1 = mat[5, ], length = 0)
  arrows(x, mat[2, ], y1 = mat[4, ], length = 0, lwd = 3)
}
ord <- order(gfra["mean", ])

# png("figures/FigSX_results_gforms.png", height = 480*10, width = 480*10, res = 72*10)
par(mai = c(1.25,0.8,0.1,0.1))
plot(c(1, 6), c(min(gfmi)+45, max(gfma)+30), type = "n", axes = FALSE, 
  ann = FALSE)
plot_ci(gfma[, ord]+5, 1:6)
plot_ci(gfmi[, ord]+45, 1:6)
plot_ci(gfra[, ord]+10, 1:6)
axis(1, at = 1:6, labels = c("Acaulescent\nrosette herbs", "Cushions", "Erect herbs", "Erect\nrosette herbs", "Prostrate\nherbs/subshrubs", "Shrubs")[ord], las = 2, cex.axis = 0.8)
axis(2, at = 0:3*5 + 40 + 5, labels = 0:3*5 + 40)
axis(2, at = 0:3*5 - 20 + 45, labels = 0:3*5 - 20)
axis(2, at = 0:4*5 + 55 + 10, labels = 0:4*5 + 55)
abline(h = 62.5, lwd = 2)
abline(h = 42.5, lwd = 2, lty = 2)
axis(2, at = c(-15+45, 47.5+5, 65+10), labels = c("Freezing resistance [°C]", 
  "Heat resistance [°C]", "TTB [°C]"), tick = FALSE, line = 2)
axis(1, labels = c("", ""), at = c(-100,100))
paste_ci <- function(mat) {
  mat <- round(mat, 2)
  paste0(mat[3, ], " [", mat[1, ], ", ", mat[5, ], "]")
}
# write.csv(rbind(paste_ci(gfra), paste_ci(gfma), paste_ci(gfmi))[, ord],
#   file = "data/analyses/gf_results.csv", row.names = FALSE)

#____________________________________________________
# localities
# png("figures/FigSX_locality.png", height = 480*10, width = 480*10, res = 72*10)
par(mai = c(1.2,0.8,0.1,0.1))
boxplot(Fr_lt50 ~ Site * Period, dat_cons, 
  col = rep(cols, 2:3),
  xlab = "", ylab = "Temperature [°C]", ylim = c(-21,25), 
  axes = FALSE, at = c(1:5,7:11,13:17,19:23))
boxplot((He_lt50-30) ~ Site * Period, dat_cons, 
  col = rep(cols, 2:3), 
  add = TRUE, at = c(1:5,7:11,13:17,19:23), axes = FALSE)
axis(2, labels = c(-50,-20,-15,-10,-5,0), at = c(-50,-20,-15,-10,-5,0))
axis(2, labels = c(35,40,45,50,55), at = c(5,10,15,20,25))
axis(1, at = c(1:5,7:11,13:17,19:23), tick = FALSE, las = 2,
  labels = rep(c("Niwot", "Brennkogel", "Antisana", "Guamani", "Yungeno"), 4))
axis(1, at = c(-100,100), labels = c("", ""))
sq <- seq(0,5, length.out = 8)
axis(2, at = sq[2:3], labels = c("",""), lwd.ticks = -1)
axis(2, at = sq[4:5], labels = c("",""), lwd.ticks = -1)
axis(2, at = sq[6:7], labels = c("",""), lwd.ticks = -1)
abline(h = mean(sq[4:5]), lty = 2)
abline(v = c(6,12,18), lty = 2)
text(c(3,9,15,21), 0, paste("Period", 1:4))

#____________________________________________________
# growth forms
# png("figures/FigSX_gforms.png", height = 480*10, width = 480*10*2, res = 72*10)
par(mai = c(1.85,0.8,0.1,0.1))
boxplot(Fr_lt50 ~ gform * Trop * Period, dat_cons, 
  col = rep(cols, each=6),
  xlab = "", ylab = "Temperature [°C]", ylim = c(-21,25), axes = FALSE, 
  at = c(1:12,14:25,27:38,40:51), las = 2)
boxplot((He_lt50-30) ~ gform * Trop * Period, dat_cons, 
  col = rep(cols, each = 6), 
  add = TRUE, at = c(1:12,14:25,27:38,40:51), axes = FALSE)
axis(2, labels = c(-50,-20,-15,-10,-5,0), at = c(-50,-20,-15,-10,-5,0))
axis(2, labels = c(35,40,45,50,55), at = c(5,10,15,20,25))
axis(1, at = c(1:12,14:25,27:38,40:51), tick = FALSE, las = 2, line = -0.5,
  labels = rep(c("Erect ros. herb", "Acaulescent ros. herb", "Erect herb", "Cushion", "Subshrub", "Shrub"), 8))
axis(1, at = c(-100,100), labels = c("", ""))
sq <- seq(0,5, length.out = 8)
axis(2, at = sq[2:3], labels = c("",""), lwd.ticks = -1)
axis(2, at = sq[4:5], labels = c("",""), lwd.ticks = -1)
axis(2, at = sq[6:7], labels = c("",""), lwd.ticks = -1)
abline(h = mean(sq[4:5]), lty = 2)
abline(v = c(13,26,39), lty = 2)
text(c(6.5,19.5,32.5,45.5), 0, paste("Period", 1:4))

#____________________________________________________
# ranges
# png("figures/FigSX_ranges.png", height = 480*10, width = 480*10, res = 72*10)
dat_cons$rng <- dat_cons$He_lt50 - dat_cons$Fr_lt50
par(mai = c(1,0.8,0.1,0.1))
boxplot(rng ~ Trop * Period, dat_cons, col = cols,
  xlab = "Period", ylab = "Temperature range [°C]", 
  at = c(2,3,5,6,8,9,11,12)/2, axes = FALSE, ylim = c(0,75), yaxs = "i")
box(bty = "l")
axis(2)
axis(1, at = 1:4*1.5-0.25, labels = 1:4, tick = FALSE)

#_