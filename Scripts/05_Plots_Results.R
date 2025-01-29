
################################################################################
################################################################################

####### Code for: The interacting demographic dynamics of the last hunter-gatherers and the early farmers (2024):
####### Cortell-Nicolau, A., Rivas, J., Crema, E. R., Shennan, S., García-Puchol, O., Kolar, J., Staniuk, R., 
####### Timpson, A.

#### Code by: Alfredo Cortell-Nicolau.

#### License: Permission is granted to use and adapt this code. Please 
#### acknowledge authorship when appropriate.

################################################################################
################                                                ################
################                   SECTION 5                    ################
################                                                ################
################################################################################

################################################################################
################                    FIGURES                     ################
################################################################################

## Load data and functions
source("Scripts/01_Functions.R")
source("Scripts/Data_for_plots.R")


################################################################################
################################################################################
#####################             MAIN FIGURES          ########################
################################################################################
################################################################################



################################################################################
#####################               FIG. 2              ########################
################################################################################

## Import libraries
library(paletteer)

## Prepare utilities
time <- 1500
colors <- colorspace::sequential_hcl(time, h = c(50,300), l = c(40,70),
                                     c = 150) ## Colour palettes
colors[1500] <- "white"
nameGn <- expression(paste(gamma["'f"])) ## expressions for the axes
nameGm <- expression(paste(gamma["hg"])) ## expressions for the axes

nameDn <- expression(paste(delta["f"])) ## expressions for the axes
nameDm <- expression(paste(delta["hg"])) ## expressions for the axes

## Indexes for looping (loading and printing)
ratio <- c(rep(0.1,3),rep(0.2,3),rep(0.4,3))
e <- rep(c(0,0.1,0.2),3)


## Values for the legend
leg_vals <- c("0-200 y","200-400 y","400-600 y","600-800 y",
              "800-1000 y","1000-1200 y","1200-1400 y","Not occurred")

## Time to peak

#png("./Figures/Fig_2.png", width = 1200, height = 850)
png("./Figures/Fig_2.png", width = 42, height = 30, unit = "cm", res = 300)
#tiff("./Figures/Fig_2.tiff", width = 1200, height = 850)

par(mfrow = c(3,3), mar = c(5,4,4,8), mai = c(0.35,0.5,0.3,0), xpd = TRUE)
for (i in 1:9){
  
  ## Load data
  dat <- readRDS(paste0("./Results/Sensitivity/Death_r_",ratio[i],"_e_",e[i],".rds"))

  ## Prepare for analytical conditions
  Dm_dat <- dat$Dm
  Gm_dat <- dat$Gm
  Gn_dat <- dat$Gn
  Dn_dat <- dat$Dn
  e_dat <- dat$e
  Km_dat <- dat$Km
  Kn_dat <- dat$Kn
  
  # Lambda values for conditions
  Li <- Gm_dat-(Dm_dat*Kn_dat)
  Lii <- Gn_dat-((Dn_dat-(e_dat*Dm_dat))*Km_dat)
  
  ## Establish parameters ranges for plots
  HG_surv <- dat[which(Li > 0 & Lii < 0),c(1,3)]
  F_surv <- dat[which(Li < 0 & Lii > 0),c(1,3)]
  Equilibrium <- dat[which(Li > 0 & Lii > 0),c(1,3)]
  One_surv <- dat[which(Li < 0 & Lii < 0),c(1,3)]
  
  ## Create convex hulls with outer points for outline polygon
  HG_hull <- chull(HG_surv$Dm,HG_surv$Dn)
  HG_hull <- c(HG_hull, HG_hull[1])
  F_hull <- chull(F_surv$Dm,F_surv$Dn)
  F_hull <- c(F_hull, F_hull[1])
  Equilibrium_hull <- chull(Equilibrium$Dm,Equilibrium$Dn)
  Equilibrium_hull <- c(Equilibrium_hull, Equilibrium_hull[1])
  One_hull <- chull(One_surv$Dm,One_surv$Dn)
  One_hull <- c(One_hull, One_hull[1])
  
  dat <- dat[,c(1,3,11)] ## 9 for Time to surpass, 10 for Time to Extinction and 11 for time to peak f
  ## Plot the whole thing
  plot(dat[,2],dat[,1],col=colors[dat[,3]], pch = 15, cex = 2.5, bty = "n", xlab = "", ylab = "")
  title(ylab = nameDm, line = 2, cex.lab = 1.5)
  mtext(nameDn, side = 1, line = 1.5, cex = 1, adj = 0.55)
  grid(lty = 2, lwd = 0.3, col = "white")
  ## Add dashed polygons
  polygon(x = HG_surv[HG_hull,2],
          y = HG_surv[HG_hull,1],
          col = adjustcolor("cyan3", alpha.f = 0.4), border = NA)
  polygon(x = F_surv[F_hull,2],
          y = F_surv[F_hull,1],
          col = adjustcolor("magenta3", alpha.f = 0.4), border = NA)
  polygon(x = Equilibrium[Equilibrium_hull,2],
          y = Equilibrium[Equilibrium_hull,1],
          col = adjustcolor("gold3", alpha.f = 0.4), border = NA)
  polygon(x = One_surv[One_hull,2],
          y = One_surv[One_hull,1],
          col = adjustcolor("lightblue3", alpha.f = 0.4), border = NA)

  # Overlay time to event
  points(dat[which(dat[,3]<1500),2], dat[which(dat[,3]<1500),1], col = colors[dat[which(dat[,3]<1500),3]], pch = 15, cex = 2.5)
  ## Add polygon squares
  lw <- 1.5
  polygon(x = HG_surv[HG_hull,2],
          y = HG_surv[HG_hull,1],
          col = NA, border = "cyan3", lwd = lw)
  polygon(x = F_surv[F_hull,2],
          y = F_surv[F_hull,1],
          col = NA, border = "magenta3", lwd = lw)
  polygon(x = Equilibrium[Equilibrium_hull,2],
          y = Equilibrium[Equilibrium_hull,1],
          col = NA, border = "gold3", lwd = lw)
  polygon(x = One_surv[One_hull,2],
          y = One_surv[One_hull,1],
          col = NA, border = "lightblue3", lwd = lw)
  
  # Add the contour (top layer)
  contour(dat[,2][1:100], dat[,1][seq(1, length(dat[,1]), 100)],
          matrix(c(dat[,3]), nrow = 100, ncol = 100), col = "black", add = TRUE)

  #if (i ==3){
  #  legend(x = 0.036, y = 0.038, legend = rep("",8), pch = rep(15,8), 
  #         col = colors[seq(1,time,length.out=8)], bty = "n",  cex = 7, y.intersp = 0.15)
  #  legend(x = 0.0435, y = 0.0335, legend = leg_vals, bty = "n",  cex = 1, y.intersp = 1.25)
  #  legend(x = 0.039, y = 0.05, legend = c("Time to","peak"), cex = 1.2, bty = "n", y.intersp = 1,
  #         text.font = 2)
  #}
  title(paste0("assimilation", " = ", e[i], " / ratio = ", ratio[i]))
}

dev.off()


################################################################################
#####################             END FIG. 2            ########################
################################################################################


################################################################################
#####################               FIG. 3              ########################
################################################################################

#png("./Figures/Fig_3.png", width = 1200, height = 850)
png("./Figures/Fig_3.png", width = 42, height = 30, unit = "cm", res = 300)
#tiff("./Figures/Fig_3.tiff", width = 1200, height = 850)

par(mfrow = c(3,3), mar = c(5,4,4,8), mai = c(0.35,0.5,0.3,0), xpd = TRUE)
for (i in 1:9){

  ## Load data
  dat <- readRDS(paste0("./Results/Sensitivity/Growth_r_",ratio[i],"_e_",e[i],".rds"))
  
  ## Prepare for analytical conditions
  Dm_dat <- dat$Dm
  Gm_dat <- dat$Gm
  Gn_dat <- dat$Gn
  Dn_dat <- dat$Dn
  e_dat <- dat$e
  Km_dat <- dat$Km
  Kn_dat <- dat$Kn
  
  # Lambda values for conditions
  Li <- Gm_dat-(Dm_dat*Kn_dat)
  Lii <- Gn_dat-((Dn_dat-(e_dat*Dm_dat))*Km_dat)
  
  ## Establish parameters ranges for plots
  HG_surv <- dat[which(Li > 0 & Lii < 0),c(2,4)]
  F_surv <- dat[which(Li < 0 & Lii > 0),c(2,4)]
  Equilibrium <- dat[which(Li > 0 & Lii > 0),c(2,4)]
  One_surv <- dat[which(Li < 0 & Lii < 0),c(2,4)]

  ## Create convex hulls with outer points for outline polygon
  HG_hull <- chull(HG_surv$Gm,HG_surv$Gn)
  HG_hull <- c(HG_hull, HG_hull[1])
  F_hull <- chull(F_surv$Gm,F_surv$Gn)
  F_hull <- c(F_hull, F_hull[1])
  Equilibrium_hull <- chull(Equilibrium$Gm,Equilibrium$Gn)
  Equilibrium_hull <- c(Equilibrium_hull, Equilibrium_hull[1])
  One_hull <- chull(One_surv$Gm,One_surv$Gn)
  One_hull <- c(One_hull, One_hull[1])
  
  dat <- dat[,c(2,4,10)] ## 9 for Time to surpass, 10 for Time to Extinction and 11 for time to peak f
  
  plot(dat[,2],dat[,1],col=colors[dat[,3]], pch = 15, cex = 2.5, bty = "n", xlab = "", ylab = "")
  title(ylab = nameGm, line = 2, cex.lab = 1.5)
  mtext(nameGn, side = 1, line = 1.5, cex = 1, adj = 0.55)
  grid(lty = 2, lwd = 0.3, col = "white")
  ## Add dashed polygons
  polygon(x = HG_surv[HG_hull,2],
          y = HG_surv[HG_hull,1],
          col = adjustcolor("cyan3", alpha.f = 0.4), border = NA)
  polygon(x = F_surv[F_hull,2],
          y = F_surv[F_hull,1],
          col = adjustcolor("magenta3", alpha.f = 0.4), border = NA)
  polygon(x = Equilibrium[Equilibrium_hull,2],
          y = Equilibrium[Equilibrium_hull,1],
          col = adjustcolor("gold3", alpha.f = 0.4), border = NA)
  polygon(x = One_surv[One_hull,2],
          y = One_surv[One_hull,1],
          col = adjustcolor("lightblue3", alpha.f = 0.4), border = NA)
  
  # Overlay time to event
  points(dat[which(dat[,3]<1500),2], dat[which(dat[,3]<1500),1], col = colors[dat[which(dat[,3]<1500),3]], pch = 15, cex = 2.5)
  ## Add polygon squares
  lw <- 1.5
  polygon(x = HG_surv[HG_hull,2],
          y = HG_surv[HG_hull,1],
          col = NA, border = "cyan3", lwd = lw)
  polygon(x = F_surv[F_hull,2],
          y = F_surv[F_hull,1],
          col = NA, border = "magenta3", lwd = lw)
  polygon(x = Equilibrium[Equilibrium_hull,2],
          y = Equilibrium[Equilibrium_hull,1],
          col = NA, border = "gold3", lwd = lw)
  polygon(x = One_surv[One_hull,2],
          y = One_surv[One_hull,1],
          col = NA, border = "lightblue3", lwd = lw)
  contour(dat[,2][1:100],dat[,1][seq(1,length(dat[,1]),100)],matrix(c(dat[,3]),nrow = 100, ncol = 100), 
          col = "black", add = TRUE)
  #if (i == 3) {
  #  legend(x = 0.0645, y = 0.019, legend = rep("",8), pch = rep(15,8), 
  #         col = colors[seq(1,time,length.out=8)], bty = "n",  cex = 7, y.intersp = 0.15)
  #  legend(x = 0.075, y = 0.0165, legend = leg_vals, bty = "n",  cex = 1, y.intersp = 1.25)
  #  legend(x = 0.0685, y = 0.02, legend = c("Time to","extinction"), cex = 1.2, bty = "n", y.intersp = 1,
  #         text.font = 2)
  #}
  
  title(paste0("assimilation", " = ", e[i], " / ratio = ", ratio[i]))
}

dev.off()



################################################################################
#####################             END FIG. 3            ########################
################################################################################



################################################################################
#####################               FIG. 4              ########################
################################################################################

## Create objects for plotting

## GENERAL

## Define x axis ticks for SPDs
plot_ticks <- list("Iberia" = data.frame("Year" = c(Neo_spd_Ib$grid$calBP[1]-1952,
                                                    Neo_spd_Ib$grid$calBP[1]-2152,
                                                    Neo_spd_Ib$grid$calBP[1]-2352,
                                                    Neo_spd_Ib$grid$calBP[1]-2552,
                                                    Neo_spd_Ib$grid$calBP[1]-2752,
                                                    Neo_spd_Ib$grid$calBP[1]-2952),
                                         "Position" =  c(2,202,402,602,802,1002)),
                   "Denmark" = data.frame("Year" = c(Neo_spd_Ja$grid$calBP[1]-1950,
                                                     Neo_spd_Ja$grid$calBP[1]-2150,
                                                     Neo_spd_Ja$grid$calBP[1]-2350,
                                                     Neo_spd_Ja$grid$calBP[1]-2550,
                                                     Neo_spd_Ja$grid$calBP[1]-2750,
                                                     Neo_spd_Ja$grid$calBP[1]-2950),
                                          "Position" =  c(0,200,400,600,800,1000)),
                   "Denmark" = data.frame("Year" = c(Neo_spd_De$grid$calBP[1]-1945,
                                                     Neo_spd_De$grid$calBP[1]-2145,
                                                     Neo_spd_De$grid$calBP[1]-2245,
                                                     Neo_spd_De$grid$calBP[1]-2345,
                                                     Neo_spd_De$grid$calBP[1]-2545,
                                                     Neo_spd_De$grid$calBP[1]-2745),
                                          "Position" =  c(-5,195,395,595,795,995)))

#png("./Figures/Fig_4.png", width = 1000, height = 700)
png("./Figures/Fig_4.png", width = 33, height = 25, unit = "cm", res = 300)
#tiff("./Figures/Fig_4.tiff", width = 1000, height = 700)

## Layout parameters
mat_plot <- matrix(c(1,1,2,3,3,4,5,5,6), nrow = 3, byrow = TRUE)
layout(mat_plot)

## Index for plot loop
Ind <- c("Ib","Ja","De")
areas <- c("Iberia","Japan","Denmark")

## Start loop
for (j in 1:length(Ind)){
  ## Create objects for plotting
  ## Assign specific object (Ib, Ja, De) to common value for loop
  Neo_spd <- get(paste0("Neo_spd_",Ind[j]))
  Meso_spd <- get(paste0("Meso_spd_",Ind[j]))
  
  fit_list <- get(paste0("fit_list_",Ind[j]))
  model_list <- get(paste0("model_list_",Ind[j]))
  
  area <- areas[j]
  
  ## Extract densities
  Dens_Neo <- sapply(fit_list, function(x) x[[3]]$grid$PrDens)
  Dens_Meso <- sapply(fit_list, function(x) x[[2]]$grid$PrDens)
  
  ## High posterior density intervals
  HPDI_Neo_low <- rep(NA,nrow(Dens_Neo))
  HPDI_Neo_high <- rep(NA,nrow(Dens_Neo))
  
  for (i in 1:nrow(Dens_Neo)){
    HPDI_Neo_low[i] <- HPDinterval(as.mcmc(Dens_Neo[i,]))[1]
    HPDI_Neo_high[i] <- HPDinterval(as.mcmc(Dens_Neo[i,]))[2]
  }
  
  HPDI_Meso_low <- rep(NA,nrow(Dens_Meso))
  HPDI_Meso_high <- rep(NA,nrow(Dens_Meso))
  
  for (i in 1:nrow(Dens_Meso)){
    HPDI_Meso_low[i] <- HPDinterval(as.mcmc(Dens_Meso[i,]))[1]
    HPDI_Meso_high[i] <- HPDinterval(as.mcmc(Dens_Meso[i,]))[2]
  }
  
  ## Extract max min vals
  Maxmin_Neo <- apply(Dens_Neo,1,maxmin)
  Maxmin_Meso <- apply(Dens_Meso,1,maxmin)
  
  ##################
  ###### Plot starts
  ##################
  
  plot(Neo_spd$grid$PrDens, type = "l", lwd = 0.1, xlab = "Years BC", ylab = "Radiocarbon density", 
       xaxt = "n", xlim = c(0,(length(HPDI_Neo_low)-1)), main = "SPD",
       ylim = c(0,max(Maxmin_Neo[2,])))
  axis(1,labels = plot_ticks[[j]]$Year, at = plot_ticks[[j]]$Position)
  legend("topleft", area[1], bty = "n", text.font = 4)
  polygon(x=c(c(0:(length(HPDI_Neo_low)-1)),c((length(HPDI_Neo_low)-1):0)), 
          y = c(HPDI_Neo_low,rev(HPDI_Neo_high)), 
          col = adjustcolor("darkslategray4", alpha.f = 0.3), border = NA)
  lines(Maxmin_Neo[1,], col = "darkslategray4", lwd = 0.7, lty = 3)
  lines(Maxmin_Neo[2,], col = "darkslategray4", lwd = 0.7, lty = 3)
  lines(Neo_spd$grid$PrDens, col = "darkslategray4", lwd = 1.5)
  
  polygon(x=c(c(0:(length(HPDI_Meso_low)-1)),c((length(HPDI_Meso_low)-1):0)), 
          y = c(HPDI_Meso_low,rev(HPDI_Meso_high)), 
          col = adjustcolor("brown4", alpha.f = 0.3), border = NA)
  lines(Maxmin_Meso[1,], col = "brown4", lwd = 0.7, lty = 3)
  lines(Maxmin_Meso[2,], col = "brown4", lwd = 0.7, lty = 3)
  lines(Meso_spd$grid$PrDens, col = "brown4", lwd = 1.5)
  
  plot(model_list[[1]][[1]][,3]/max(model_list[[1]][[1]][,3]), ylim = c(0,1),
       col = adjustcolor("darkslategray4", alpha.f = 0.05), type = "l",
       ylab = "Population density", xlab = "Time", main = "Population")
  for (i in 2:500){
    lines(model_list[[i]][[1]][,3]/max(model_list[[i]][[1]][,3]), 
          col = adjustcolor("darkslategray4", alpha.f = 0.05), type = "l")
  }
  for (i in 1:500){
    lines(model_list[[i]][[1]][,2]/max(model_list[[i]][[1]][,3]), 
          col = adjustcolor("brown4", alpha.f = 0.05), type = "l")
  }
}

dev.off()

################################################################################
#####################             END FIG. 4            ########################
################################################################################


################################################################################
#####################               FIG. 5              ########################
################################################################################

library(vioplot)

## Prepare utilities
params_Ib <- data.frame("Dm" = sapply(model_list_Ib, function(x) x[[2]]$Dm),
                        "Gm" = sapply(model_list_Ib, function(x) x[[2]]$Gm),
                        "Dn" = sapply(model_list_Ib, function(x) x[[2]]$Dn),
                        "Gn" = sapply(model_list_Ib, function(x) x[[2]]$Gn),
                        "e" = sapply(model_list_Ib, function(x) x[[2]]$e),
                        "m" = sapply(model_list_Ib, function(x) x[[2]]$m))
params_Ib$Gn <- params_Ib$Gn+params_Ib$m

params_Ja <- data.frame("Dm" = sapply(model_list_Ja, function(x) x[[2]]$Dm),
                        "Gm" = sapply(model_list_Ja, function(x) x[[2]]$Gm),
                        "Dn" = sapply(model_list_Ja, function(x) x[[2]]$Dn),
                        "Gn" = sapply(model_list_Ja, function(x) x[[2]]$Gn),
                        "e" = sapply(model_list_Ja, function(x) x[[2]]$e),
                        "m" = sapply(model_list_Ja, function(x) x[[2]]$m))
params_Ja$Gn <- params_Ja$Gn+params_Ja$m

params_De <- data.frame("Dm" = sapply(model_list_De, function(x) x[[2]]$Dm),
                        "Gm" = sapply(model_list_De, function(x) x[[2]]$Gm),
                        "Dn" = sapply(model_list_De, function(x) x[[2]]$Dn),
                        "Gn" = sapply(model_list_De, function(x) x[[2]]$Gn),
                        "e" = sapply(model_list_De, function(x) x[[2]]$e),
                        "m" = sapply(model_list_De, function(x) x[[2]]$m))
params_De$Gn <- params_De$Gn+params_De$m

# Index to select parameters
sel_par <- c("Dm","Gm","Dn","Gn","e")

## To assign names with math letters
nameDm <- expression(paste(delta["hg"]))
nameGm <- expression(paste(gamma["hg"]))
nameDn <- expression(paste(delta["f"]))
nameGn <- expression(paste(gamma["'f"]))
namee <- expression(eta)

#png("./Figures/Fig_5.png", width = 1000, height = 700)
png("./Figures/Fig_5.png", width = 33, height = 25, unit = "cm", res = 300)
#tiff("./Figures/Fig_5.tiff", width = 1000, height = 700)

par(mfrow = c(2,3))

## Prepare density for legend
normdens <- density(rnorm(10000,32,0.5), bw = 4)

plot(x=normdens$x,y=normdens$y/max(normdens$y)+5, type = "l", ylim = c(1,7.5), xlim = c(0,50),
     col = "white", axes = FALSE, ylab = "", xlab = "")
polygon(normdens$x,normdens$y/max(normdens$y)+5, col = adjustcolor("mediumorchid4", alpha.f = 0.6), border = "gray")
polygon(normdens$x,normdens$y/max(normdens$y)+3, col = adjustcolor("cadetblue3", alpha.f = 0.6), border = "gray")
polygon(normdens$x,normdens$y/max(normdens$y)+1, col = adjustcolor("sienna4", alpha.f = 0.6), border = "gray")

mtext("Denmark", 2, at = 5.5, las = 2, line = -4.5, col = "gray25", cex = 1.6, font = 2)
mtext("Iberia", 2, at = 3.5, las = 2, line = -4.5, col = "gray25", cex = 1.6, font = 2)
mtext("Japan", 2, at = 1.5, las = 2, line = -4.5, col = "gray25", cex = 1.6, font = 2)

for (i in 1:5){
  parDe <- params_De[,i]
  parIb <- params_Ib[,i]
  parJa <- params_Ja[,i]
  
  ## To automatically create ylimit depending on the area
  ylobj <- c(parDe,parIb,parJa)
  yl <- max(ylobj) + (max(ylobj)-quantile(max(ylobj))[2])
  
  # Common parameters
  peach <- 21
  alph <- 0.4
  cexx <- 0.6
  
  set.seed(1)
  
  ## Create structure and background
  plot(x = c(1,7), y = c(0,0), type = "l", col = "gray", lwd = 0, axes = FALSE, ylab = "Density", xlab = "",
       ylim = c(0,yl))
  line_bg <- axis(2, tck = -0.015)
  abline(h=line_bg,col = "gray",lwd = 0.5)
  title(get(paste0("name",sel_par[i])), cex.main = 2)
  
  ## Plot distributions
  ## Denmark
  points(x = rnorm(length(parDe),2,0.2), parDe, col = adjustcolor("black", alpha.f = alph), 
         pch = peach, cex = cexx, bg = adjustcolor("mediumorchid4", alpha.f = alph))
  vioplot(parDe, col = adjustcolor("mediumorchid4", alpha.f = 0.6), border = "gray", add = TRUE, drawRect = FALSE,
          at = 2, pchMed = NA, colMed = NA)
  rect(xleft = 1.95, ybottom = HPDinterval(as.mcmc(parDe))[,1], xright = 2.05, ytop = HPDinterval(as.mcmc(parDe))[,2], 
       col = "mediumorchid4", border = "gray")
  segments(x1 = 1, y1 = median(parDe),x0 = 3, y0 = median(parDe), col = "red", lty = 2, lwd = 1.2)
  
  ## Iberia
  points(x = rnorm(length(parIb),4,0.2), parIb, col = adjustcolor("black", alpha.f = alph), pch = peach, cex = cexx, 
         bg = adjustcolor("cadetblue3", alpha.f = alph))
  vioplot(parIb, col = adjustcolor("cadetblue3", alpha.f = 0.6), border = "gray", add = TRUE, drawRect = FALSE,
          at = 4, pchMed = NA, colMed = NA)
  rect(xleft = 3.95, ybottom = HPDinterval(as.mcmc(parIb))[,1], xright = 4.05, ytop = HPDinterval(as.mcmc(parIb))[,2], 
       col = "cadetblue4", border = "gray")
  segments(x1 = 3, y1 = median(parIb),x0 = 5, y0 = median(parIb), col = "red", lty = 2, lwd = 1.2)
  
  
  ## Japan
  points(x = rnorm(length(parJa),6,0.2), parJa, col = adjustcolor("black", alpha.f = alph), pch = peach, cex = cexx, 
         bg = adjustcolor("sienna4", alpha.f = alph))
  vioplot(parJa, col = adjustcolor("sienna4", alpha.f = 0.6), border = "gray", add = TRUE, drawRect = FALSE,
          at = 6, pchMed = NA, colMed = NA)
  rect(xleft = 5.95, ybottom = HPDinterval(as.mcmc(parJa))[,1], xright = 6.05, ytop = HPDinterval(as.mcmc(parJa))[,2], 
       col = "sienna4", border = "gray")
  segments(x1 = 5, y1 = median(parJa),x0 = 7, y0 = median(parJa), col = "red", lty = 2, lwd = 1.2)
  
}

dev.off()


################################################################################
#####################             END FIG. 5            ########################
################################################################################


################################################################################
#####################               FIG. 6              ########################
################################################################################

### This includes the ratio comparison

### COMPARE RATIOS

ext_ratios <- function(x){
  iPM <- x[[1]][,2][1]
  iF <- x[[1]][,3][1]
  
  return(iF/iPM)
}

## Extract ratios
In_ratio_De <- sapply(model_list_De, ext_ratios)
In_ratio_Ja <- sapply(model_list_Ja, ext_ratios)
In_ratio_Ib <- sapply(model_list_Ib, ext_ratios)

## Means, medians and HPDIs
## Denmark
mean_ratio_De <- mean(In_ratio_De)
median_ratio_De <- median(In_ratio_De)
HPDI_ratio_De <- HPDinterval(as.mcmc(In_ratio_De))

## Japan
mean_ratio_Ja <- mean(In_ratio_Ja)
median_ratio_Ja <- median(In_ratio_Ja)
HPDI_ratio_Ja <- HPDinterval(as.mcmc(In_ratio_Ja))

## Iberia
mean_ratio_Ib <- mean(In_ratio_Ib)
median_ratio_Ib <- median(In_ratio_Ib)
HPDI_ratio_Ib <- HPDinterval(as.mcmc(In_ratio_Ib))

ratio_res <- data.frame("Denmark" = c(mean_ratio_De,median_ratio_De,HPDI_ratio_De),
                        "Japan" = c(mean_ratio_Ja,median_ratio_Ja,HPDI_ratio_Ja),
                        "Iberia" = c(mean_ratio_Ib,median_ratio_Ib,HPDI_ratio_Ib))
rownames(ratio_res) <- c("mean", "median", "low.HPDI", "high.HPDI")
#View(ratio_res)


### TIMES TO DIFFERENT EVENTS

## Time to K not compared because the case-studies were designed to stop when F reaches K

## Time to surass
time_to_surpass <- function(x){
  obj <- x[[1]]
  t_t_s <- obj[obj[,3]>obj[,2],][1,1]
  return(t_t_s)
}
t_t_s_De <- sapply(model_list_De, time_to_surpass)
t_t_s_Ja <- sapply(model_list_Ja, time_to_surpass)
t_t_s_Ib <- sapply(model_list_Ib, time_to_surpass)

## Time to disappear
time_to_dis <- function(x){
  obj <- x[[1]]
  t_t_d <- obj[obj[,2]<= obj[1,3]/20,][1,1] ## Extinction when hg population is lower than 1/20 of initial farmer population
  #ocurrence <- table(is.na(t_t_d))[2]
  return(t_t_d)
}

t_t_d_De <- sapply(model_list_De, time_to_dis)
t_t_d_Ja <- sapply(model_list_Ja, time_to_dis)
t_t_d_Ib <- sapply(model_list_Ib, time_to_dis)

## Extract HPDIs
t_t_s_De_HPDI_low <- HPDinterval(as.mcmc(t_t_s_De))[1]
t_t_s_De_HPDI_high <- HPDinterval(as.mcmc(t_t_s_De))[2]
t_t_s_Ja_HPDI_low <- HPDinterval(as.mcmc(t_t_s_Ja))[1]
t_t_s_Ja_HPDI_high <- HPDinterval(as.mcmc(t_t_s_Ja))[2]
t_t_s_Ib_HPDI_low <- HPDinterval(as.mcmc(t_t_s_Ib))[1]
t_t_s_Ib_HPDI_high <- HPDinterval(as.mcmc(t_t_s_Ib))[2]

t_t_d_De_HPDI_low <- HPDinterval(as.mcmc(t_t_d_De))[1]
t_t_d_De_HPDI_high <- HPDinterval(as.mcmc(t_t_d_De))[2]
t_t_d_Ja_HPDI_low <- HPDinterval(as.mcmc(t_t_d_Ja))[1]
t_t_d_Ja_HPDI_high <- HPDinterval(as.mcmc(t_t_d_Ja))[2]
t_t_d_Ib_HPDI_low <- HPDinterval(as.mcmc(t_t_d_Ib))[1]
t_t_d_Ib_HPDI_high <- HPDinterval(as.mcmc(t_t_d_Ib))[2]

## Prepare densities
dens_t_t_s_De <- density(t_t_s_De)
dens_t_t_s_De_df <- data.frame("x" = dens_t_t_s_De$x, "y" = dens_t_t_s_De$y)
dens_t_t_s_De_df_HPDI <- dens_t_t_s_De_df[dens_t_t_s_De_df$x >= t_t_s_De_HPDI_low & dens_t_t_s_De_df$x <= t_t_s_De_HPDI_high,]

dens_t_t_s_Ja <- density(t_t_s_Ja)
dens_t_t_s_Ja_df <- data.frame("x" = dens_t_t_s_Ja$x, "y" = dens_t_t_s_Ja$y)
dens_t_t_s_Ja_df_HPDI <- dens_t_t_s_Ja_df[dens_t_t_s_Ja_df$x >= t_t_s_Ja_HPDI_low & dens_t_t_s_Ja_df$x <= t_t_s_Ja_HPDI_high,]

dens_t_t_s_Ib <- density(t_t_s_Ib)
dens_t_t_s_Ib_df <- data.frame("x" = dens_t_t_s_Ib$x, "y" = dens_t_t_s_Ib$y)
dens_t_t_s_Ib_df_HPDI <- dens_t_t_s_Ib_df[dens_t_t_s_Ib_df$x >= t_t_s_Ib_HPDI_low & dens_t_t_s_Ib_df$x <= t_t_s_Ib_HPDI_high,]

dens_t_t_d_De <- density(na.omit(t_t_d_De))
dens_t_t_d_De_df <- data.frame("x" = dens_t_t_d_De$x, "y" = dens_t_t_d_De$y)
dens_t_t_d_De_df_HPDI <- dens_t_t_d_De_df[dens_t_t_d_De_df$x >= t_t_d_De_HPDI_low & dens_t_t_d_De_df$x <= t_t_d_De_HPDI_high,]

dens_t_t_d_Ja <- density(na.omit(t_t_d_Ja))
dens_t_t_d_Ja_df <- data.frame("x" = dens_t_t_d_Ja$x, "y" = dens_t_t_d_Ja$y)
dens_t_t_d_Ja_df_HPDI <- dens_t_t_d_Ja_df[dens_t_t_d_Ja_df$x >= t_t_d_Ja_HPDI_low & dens_t_t_d_Ja_df$x <= t_t_d_Ja_HPDI_high,]

dens_t_t_d_Ib <- density(na.omit(t_t_d_Ib))
dens_t_t_d_Ib_df <- data.frame("x" = dens_t_t_d_Ib$x, "y" = dens_t_t_d_Ib$y)
dens_t_t_d_Ib_df_HPDI <- dens_t_t_d_Ib_df[dens_t_t_d_Ib_df$x >= t_t_d_Ib_HPDI_low & dens_t_t_d_Ib_df$x <= t_t_d_Ib_HPDI_high,]


## Plots
#png("./Figures/Fig_6.png", width = 1200, height = 850)
png("./Figures/Fig_6.png", width = 40, height = 30, unit = "cm", res = 300)
#tiff("./Figures/Fig_6.tiff", width = 1200, height = 850)

par(mfrow = c(3,2))

## Time to surpass
plot(dens_t_t_s_De, main = "Time to surpass", col = "lightblue4", ylab = "", xlab = "Years", axes = FALSE, xlim = c(0,800))
polygon(dens_t_t_s_De_df$x, dens_t_t_s_De_df$y, col = adjustcolor("lightblue", alpha.f = 0.4), border = NA)
polygon(c(dens_t_t_s_De_df_HPDI$x,rev(dens_t_t_s_De_df_HPDI$x)),
        c(rep(0.00001,length(dens_t_t_s_De_df_HPDI$x)),rev(dens_t_t_s_De_df_HPDI$y)), 
        col = adjustcolor("lightblue4", alpha.f = 1), border = NA)
abline(v = median(t_t_s_De), lty = 2, col = "white")
axis(1)
mtext("Denmark",2, cex = 1.2, font = 2, line = 1)

plot(dens_t_t_d_De, main = "Time to hg disappearance", col = "lightblue4", ylab = "", xlab = "Years", axes = FALSE, xlim = c(0,800))
polygon(dens_t_t_d_De_df$x, dens_t_t_d_De_df$y, col = adjustcolor("lightblue", alpha.f = 0.4), border = NA)
polygon(c(dens_t_t_d_De_df_HPDI$x,rev(dens_t_t_d_De_df_HPDI$x)),
        c(rep(0.00001,length(dens_t_t_d_De_df_HPDI$x)),rev(dens_t_t_d_De_df_HPDI$y)), 
        col = adjustcolor("lightblue4", alpha.f = 1), border = NA)
abline(v = median(na.omit(t_t_d_De)), lty = 2, col = "white")
axis(1)

plot(dens_t_t_s_Ja, col = "lightblue4", ylab = "", xlab = "Years", main = "", axes = FALSE, xlim = c(0,800))
polygon(dens_t_t_s_Ja_df$x, dens_t_t_s_Ja_df$y, col = adjustcolor("lightblue", alpha.f = 0.4), border = NA)
polygon(c(dens_t_t_s_Ja_df_HPDI$x,rev(dens_t_t_s_Ja_df_HPDI$x)),
        c(rep(0.00001,length(dens_t_t_s_Ja_df_HPDI$x)),rev(dens_t_t_s_Ja_df_HPDI$y)), 
        col = adjustcolor("lightblue4", alpha.f = 1), border = NA)
abline(v = median(t_t_s_Ja), lty = 2, col = "white")
axis(1)
mtext("Japan",2, cex = 1.2, font = 2, line = 1)

plot(dens_t_t_d_Ja, col = "lightblue4", ylab = "", xlab = "Years", main = "", axes = FALSE, xlim = c(0,800))
polygon(dens_t_t_d_Ja_df$x, dens_t_t_d_Ja_df$y, col = adjustcolor("lightblue", alpha.f = 0.4), border = NA)
polygon(c(dens_t_t_d_Ja_df_HPDI$x,rev(dens_t_t_d_Ja_df_HPDI$x)),
        c(rep(0.00001,length(dens_t_t_d_Ja_df_HPDI$x)),rev(dens_t_t_d_Ja_df_HPDI$y)), 
        col = adjustcolor("lightblue4", alpha.f = 1), border = NA)
abline(v = median(t_t_d_Ja), lty = 2, col = "white")
axis(1)

plot(dens_t_t_s_Ib, main = "", col = "lightblue4", ylab = "", xlab = "Years", axes = FALSE, xlim = c(0,800))
polygon(dens_t_t_s_Ib_df$x, dens_t_t_s_Ib_df$y, col = adjustcolor("lightblue", alpha.f = 0.4), border = NA)
polygon(c(dens_t_t_s_Ib_df_HPDI$x,rev(dens_t_t_s_Ib_df_HPDI$x)),
        c(rep(0.00001,length(dens_t_t_s_Ib_df_HPDI$x)),rev(dens_t_t_s_Ib_df_HPDI$y)), 
        col = adjustcolor("lightblue4", alpha.f = 1), border = NA)
abline(v = median(t_t_s_Ib), lty = 2, col = "white")
axis(1)
mtext("Iberia",2, cex = 1.2, font = 2, line = 1)


plot(dens_t_t_d_Ib, main = "", col = "lightblue4", ylab = "", xlab = "Years", axes = FALSE, xlim = c(0,800))
polygon(dens_t_t_d_Ib_df$x, dens_t_t_d_Ib_df$y, col = adjustcolor("lightblue", alpha.f = 0.4), border = NA)
polygon(c(dens_t_t_d_Ib_df_HPDI$x,rev(dens_t_t_d_Ib_df_HPDI$x)),
        c(rep(0.00001,length(dens_t_t_d_Ib_df_HPDI$x)),rev(dens_t_t_d_Ib_df_HPDI$y)), 
        col = adjustcolor("lightblue4", alpha.f = 1), border = NA)
abline(v = median(na.omit(t_t_d_Ib)), lty = 2, col = "white")
axis(1)

dev.off()


################################################################################
#####################             END FIG. 6            ########################
################################################################################



################################################################################
################################################################################
#####################        FIGURES SUPLEMENTARY       ########################
################################################################################
################################################################################



################################################################################
#####################              FIG. S1 1            ########################
################################################################################
state <- c(M = 60, N = 6)
methds1 <- readRDS("./Data/pars_lilower0_lihigher0.rds")
out <- data.frame(LotVolme(methds1, times = c(1:600)))

#png("./Figures/Fig_S1_1.png", width = 1000, height = 700)
png("./Figures/Fig_S1_1.png", width = 33, height = 25, unit = "cm", res = 300)
#tiff("./Figures/Fig_S1_1.tiff", width = 1000, height = 700)

plot(out[,2], type = "l", ylim = c(0,max(out[,3])), col = "black", lwd = 1.5,
     ylab = "Population", xlab = "time", main = expression(paste(lambda["1"]," < 0 and ", lambda["2"], " > 0")))
lines(out[,3], col = "red", lwd = 1.5)

dev.off()

################################################################################
#####################            END FIG. S1 1          ########################
################################################################################



################################################################################
#####################              FIG. S1 2            ########################
################################################################################

methds2 <- readRDS("./Data/pars_lilower0_lilower0.rds")
out <- data.frame(LotVolme(methds2, times = c(1:600)))

#png("./Figures/Fig_S1_2.png", width = 1000, height = 700)
png("./Figures/Fig_S1_2.png", width = 33, height = 25, unit = "cm", res = 300)
#tiff("./Figures/Fig_S1_2.tiff", width = 1000, height = 700)

plot(out[,2], type = "l", ylim = c(0,max(out[,2])), col = "black", lwd = 1.5,
     ylab = "Population", xlab = "time", main = expression(paste(lambda["1"]," < 0 and ", lambda["2"], " < 0")))
lines(out[,3], col = "red", lwd = 1.5)

dev.off()

################################################################################
#####################            END FIG. S1 1          ########################
################################################################################


################################################################################
#####################              FIG. S2 1            ########################
################################################################################

## Total number of dates
Tot_dates <- Neo_spd_sim$metadata$ndates + Meso_spd_sim$metadata$ndates

## Extract densities
Dens_Neo <- sapply(fit_list_sim, function(x) x[[3]]$grid$PrDens)
Dens_Meso <- sapply(fit_list_sim, function(x) x[[2]]$grid$PrDens)

## High posterior density intervals
HPDI_Neo_low <- rep(NA,nrow(Dens_Neo))
HPDI_Neo_high <- rep(NA,nrow(Dens_Neo))

for (i in 1:nrow(Dens_Neo)){
  HPDI_Neo_low[i] <- HPDinterval(as.mcmc(Dens_Neo[i,]))[1]
  HPDI_Neo_high[i] <- HPDinterval(as.mcmc(Dens_Neo[i,]))[2]
}

HPDI_Meso_low <- rep(NA,nrow(Dens_Meso))
HPDI_Meso_high <- rep(NA,nrow(Dens_Meso))

for (i in 1:nrow(Dens_Meso)){
  HPDI_Meso_low[i] <- HPDinterval(as.mcmc(Dens_Meso[i,]))[1]
  HPDI_Meso_high[i] <- HPDinterval(as.mcmc(Dens_Meso[i,]))[2]
}

## Extract max min vals
Maxmin_Neo <- apply(Dens_Neo,1,maxmin)
Maxmin_Meso <- apply(Dens_Meso,1,maxmin)

## Plot figure

#png("./Figures/Fig_S2_1.png", width = 1000, height = 700)
png("./Figures/Fig_S2_1.png", width = 33, height = 25, unit = "cm", res = 300)
#tiff("./Figures/Fig_S2_1.tiff", width = 1000, height = 700)

mat_layout <- matrix(c(rep(1,6),seq(2,7)),byrow = TRUE,nrow = 4)
layout(mat_layout)

plot(Neo_spd_sim$grid$PrDens, type = "l", lwd = 0.1, xlab = "Years BC", ylab = "Density", 
     xaxt = "n", xlim = c(0,length(HPDI_Neo_low)-1), main = "Simulated",
     ylim = c(0,max(Maxmin_Neo[2,])))
axis(1,labels = c(Neo_spd_sim$grid$calBP[1]-1950,
                  Neo_spd_sim$grid$calBP[1]-2150,
                  Neo_spd_sim$grid$calBP[1]-2350,
                  Neo_spd_sim$grid$calBP[1]-2550,
                  Neo_spd_sim$grid$calBP[1]-2750,
                  Neo_spd_sim$grid$calBP[1]-2950), at = c(0,200,400,600,800,1000))
axis(3,labels = paste0("ndates = ", Tot_dates), at = 50, tick = F, font = 4)
polygon(x=c(c(0:(length(HPDI_Neo_low)-1)),c((length(HPDI_Neo_low)-1):0)), 
        y = c(HPDI_Neo_low,rev(HPDI_Neo_high)), 
        col = adjustcolor("steelblue4", alpha.f = 0.3), border = NA)
lines(Maxmin_Neo[1,], col = "steelblue4", lwd = 0.7, lty = 3)
lines(Maxmin_Neo[2,], col = "steelblue4", lwd = 0.7, lty = 3)
lines(Neo_spd_sim$grid$PrDens, col = "steelblue4", lwd = 1.5)

polygon(x=c(c(0:(length(HPDI_Meso_low)-1)),c((length(HPDI_Meso_low)-1):0)), 
        y = c(HPDI_Meso_low,rev(HPDI_Meso_high)), 
        col = adjustcolor("lightgoldenrod4", alpha.f = 0.3), border = NA)
lines(Maxmin_Meso[1,], col = "lightgoldenrod4", lwd = 0.7, lty = 3)
lines(Maxmin_Meso[2,], col = "lightgoldenrod4", lwd = 0.7, lty = 3)
lines(Meso_spd_sim$grid$PrDens, col = "lightgoldenrod4", lwd = 1.5)

################################################################################
################################################################################
################################################################################

## PARAMETERS
# Simulated parameters
params <- data.frame("Dm" = sapply(model_list_sim, function(x) x[[2]]$Dm),
                     "Gm" = sapply(model_list_sim, function(x) x[[2]]$Gm),
                     "Dn" = sapply(model_list_sim, function(x) x[[2]]$Dn),
                     "Gn" = sapply(model_list_sim, function(x) x[[2]]$Gn),
                     "e" = sapply(model_list_sim, function(x) x[[2]]$e),
                     "m" = sapply(model_list_sim, function(x) x[[2]]$m))
params$`'n` <- params$Gn+params$m

## Target parameters
Tact_pars <- c("Dm" = 0.002,
               "Gm" = 0.002,
               "Dn" = 0.0015,
               "Gn" = 0.025,
               "e" = 0.2,
               "m" = 0.01,
               "'n" = 0.035)

# Index to select parameters
sel_par <- c("Dm","Gm","Dn","Gn","e")

## To assign names with math letters
nameDm <- expression(paste(delta["hg"]))
nameGm <- expression(paste(gamma["hg"]))
nameDn <- expression(paste(delta["f"]))
nameGn <- expression(paste(gamma["'f"]))
namee <- expression(eta)
#par(mfrow = c(2,3))

for (i in 1:5){
  Dens <- density(params[,sel_par[i]], adjust = 2)
  
  ## Compute credible intervals to show in plot
  HPDI_95 <- HPDinterval(as.mcmc(params[,sel_par[i]]))
  HPDI_95_range <- Dens$x > HPDI_95[1] & Dens$x < HPDI_95[2] 
  
  HPDI_80 <- HPDinterval(as.mcmc(params[,sel_par[i]]), 0.8)
  HPDI_80_range <- Dens$x > HPDI_80[1] & Dens$x < HPDI_80[2] 
  
  plot(Dens, lwd = 0.1, main = get(paste0("name",sel_par[i])), xlab = "value", 
       ylim = c(0,max(Dens$y)+5))
  polygon(c(HPDI_95[1], Dens$x[HPDI_95_range], HPDI_95[2]),                
          c(0, Dens$y[HPDI_95_range],0),
          col = adjustcolor("steelblue4", alpha.f = 0.5), border = "black", lty = 3,
          lwd = 0.5)
  polygon(c(HPDI_80[1], Dens$x[HPDI_80_range], HPDI_80[2]),                
          c(0, Dens$y[HPDI_80_range],0),
          col = adjustcolor("steelblue", alpha.f = 0.5), border = "black", lty = 3,
          lwd = 0.5)
  abline(v=Tact_pars[sel_par[i]], lty = 3, col = "black", lwd = 0.8)
  abline(v=mean(params[,sel_par[i]]), lty = 3, col = "red", lwd = 0.8)
}

## To report if needed
means_params <- c("Dm" = mean(params$Dm),
                  "Gm" = mean(params$Gm),
                  "Dn" = mean(params$Dn),
                  "Gn" = mean(params$Gn),
                  "e" = mean(params$e),
                  "m" = mean(params$m),
                  "'n" = mean(params$`'n`))

vars_params <- c("Dm" = var(params$Dm),
                 "Gm" = var(params$Gm),
                 "Dn" = var(params$Dn),
                 "Gn" = var(params$Gn),
                 "e" = var(params$e),
                 "m" = var(params$m),
                 "'n" = var(params$`'n`))


dev.off()

################################################################################
#####################           END FIG. S2 1           ########################
################################################################################


################################################################################
#####################              FIG. S3 3            ########################
################################################################################

## Prepare utilities
time <- 1500
colors <- colorspace::sequential_hcl(time, h = c(50,300), l = c(40,70),
                                     c = 150) ## Colour palettes
colors[1500] <- "white"
nameGn <- expression(paste(gamma["'f"])) ## expressions for the axes
nameGm <- expression(paste(gamma["hg"])) ## expressions for the axes

nameDn <- expression(paste(delta["f"])) ## expressions for the axes
nameDm <- expression(paste(delta["hg"])) ## expressions for the axes

## Indexes for looping (loading and printing)
ratio <- c(rep(0.1,3),rep(0.2,3),rep(0.4,3))
e <- rep(c(0,0.1,0.2),3)


## Values for the legend
leg_vals <- c("0-200 y","200-400 y","400-600 y","600-800 y",
              "800-1000 y","1000-1200 y","1200-1400 y","Not occurred")

## Time to surpass

#png("./Figures/Fig_S3_3.png", width = 1200, height = 850)
png("./Figures/Fig_S3_3.png", width = 41, height = 30, unit = "cm", res = 300)
#tiff("./Figures/Fig_S3_3.tiff", width = 1200, height = 850)

par(mfrow = c(3,3), mar = c(5,4,4,8), mai = c(0.35,0.5,0.3,0), xpd = TRUE)
for (i in 1:9){
  
  ## Load data
  dat <- readRDS(paste0("./Results/Sensitivity/Death_r_",ratio[i],"_e_",e[i],".rds"))
  dat$time.to.overcome[which(dat$time.to.overcome == 1501)] <- 1500
  
  
  ## Prepare for analytical conditions
  Dm_dat <- dat$Dm
  Gm_dat <- dat$Gm
  Gn_dat <- dat$Gn
  Dn_dat <- dat$Dn
  e_dat <- dat$e
  Km_dat <- dat$Km
  Kn_dat <- dat$Kn
  
  # Lambda values for conditions
  Li <- Gm_dat-(Dm_dat*Kn_dat)
  Lii <- Gn_dat-((Dn_dat-(e_dat*Dm_dat))*Km_dat)
  
  ## Establish parameters ranges for plots
  HG_surv <- dat[which(Li > 0 & Lii < 0),c(1,3)]
  F_surv <- dat[which(Li < 0 & Lii > 0),c(1,3)]
  Equilibrium <- dat[which(Li > 0 & Lii > 0),c(1,3)]
  One_surv <- dat[which(Li < 0 & Lii < 0),c(1,3)]
  
  ## Create convex hulls with outer points for outline polygon
  HG_hull <- chull(HG_surv$Dm,HG_surv$Dn)
  HG_hull <- c(HG_hull, HG_hull[1])
  F_hull <- chull(F_surv$Dm,F_surv$Dn)
  F_hull <- c(F_hull, F_hull[1])
  Equilibrium_hull <- chull(Equilibrium$Dm,Equilibrium$Dn)
  Equilibrium_hull <- c(Equilibrium_hull, Equilibrium_hull[1])
  One_hull <- chull(One_surv$Dm,One_surv$Dn)
  One_hull <- c(One_hull, One_hull[1])
  
  dat <- dat[,c(1,3,9)] ## 9 for Time to surpass, 10 for Time to Extinction and 11 for time to peak f
  
  plot(dat[,2],dat[,1],col=colors[dat[,3]], pch = 15, cex = 2.5, bty = "n", xlab = "", ylab = "")
  title(ylab = nameDm, line = 2, cex.lab = 1.5)
  mtext(nameDn, side = 1, line = 1.5, cex = 1, adj = 0.55)
  grid(lty = 2, lwd = 0.3, col = "white")
  ## Add dashed polygons
  polygon(x = HG_surv[HG_hull,2],
          y = HG_surv[HG_hull,1],
          col = adjustcolor("cyan3", alpha.f = 0.4), border = NA)
  polygon(x = F_surv[F_hull,2],
          y = F_surv[F_hull,1],
          col = adjustcolor("magenta3", alpha.f = 0.4), border = NA)
  polygon(x = Equilibrium[Equilibrium_hull,2],
          y = Equilibrium[Equilibrium_hull,1],
          col = adjustcolor("gold3", alpha.f = 0.4), border = NA)
  polygon(x = One_surv[One_hull,2],
          y = One_surv[One_hull,1],
          col = adjustcolor("lightblue3", alpha.f = 0.4), border = NA)
  
  # Overlay time to event
  points(dat[which(dat[,3]<1500),2], dat[which(dat[,3]<1500),1], col = colors[dat[which(dat[,3]<1500),3]], pch = 15, cex = 2.5)
  ## Add polygon squares
  lw <- 1.5
  polygon(x = HG_surv[HG_hull,2],
          y = HG_surv[HG_hull,1],
          col = NA, border = "cyan3", lwd = lw)
  polygon(x = F_surv[F_hull,2],
          y = F_surv[F_hull,1],
          col = NA, border = "magenta3", lwd = lw)
  polygon(x = Equilibrium[Equilibrium_hull,2],
          y = Equilibrium[Equilibrium_hull,1],
          col = NA, border = "gold3", lwd = lw)
  polygon(x = One_surv[One_hull,2],
          y = One_surv[One_hull,1],
          col = NA, border = "lightblue3", lwd = lw)
  contour(dat[,2][1:100],dat[,1][seq(1,length(dat[,1]),100)],matrix(c(dat[,3]),nrow = 100, ncol = 100), 
          col = "black", add = TRUE)
  #if (i == 3){
  #  legend(x = 0.036, y = 0.038, legend = rep("",8), pch = rep(15,8), 
  #         col = colors[seq(1,time,length.out=8)], bty = "n",  cex = 7, y.intersp = 0.15)
  #  legend(x = 0.0435, y = 0.0335, legend = leg_vals, bty = "n",  cex = 1, y.intersp = 1.25)
  #  legend(x = 0.039, y = 0.05, legend = c("Time to","surpass"), cex = 1.2, bty = "n", y.intersp = 1,
  #         text.font = 2)
  #}
  title(paste0("assimilation", " = ", e[i], " / ratio = ", ratio[i]))
}

dev.off()


################################################################################
#####################           END FIG. S3 3           ########################
################################################################################



################################################################################
#####################              FIG. S3 4            ########################
################################################################################

## Time to extinction

#png("./Figures/Fig_S3_4.png", width = 1200, height = 850)
png("./Figures/Fig_S3_4.png", width = 41, height = 30, unit = "cm", res = 300)
#tiff("./Figures/Fig_S3_4.tiff", width = 1200, height = 850)

par(mfrow = c(3,3), mar = c(5,4,4,8), mai = c(0.35,0.5,0.3,0), xpd = TRUE)
for (i in 1:9){
  
  ## Load data
  dat <- readRDS(paste0("./Results/Sensitivity/Death_r_",ratio[i],"_e_",e[i],".rds"))
  
  ## Prepare for analytical conditions
  Dm_dat <- dat$Dm
  Gm_dat <- dat$Gm
  Gn_dat <- dat$Gn
  Dn_dat <- dat$Dn
  e_dat <- dat$e
  Km_dat <- dat$Km
  Kn_dat <- dat$Kn
  
  # Lambda values for conditions
  Li <- Gm_dat-(Dm_dat*Kn_dat)
  Lii <- Gn_dat-((Dn_dat-(e_dat*Dm_dat))*Km_dat)
  
  ## Establish parameters ranges for plots
  HG_surv <- dat[which(Li > 0 & Lii < 0),c(1,3)]
  F_surv <- dat[which(Li < 0 & Lii > 0),c(1,3)]
  Equilibrium <- dat[which(Li > 0 & Lii > 0),c(1,3)]
  One_surv <- dat[which(Li < 0 & Lii < 0),c(1,3)]
  
  ## Create convex hulls with outer points for outline polygon
  HG_hull <- chull(HG_surv$Dm,HG_surv$Dn)
  HG_hull <- c(HG_hull, HG_hull[1])
  F_hull <- chull(F_surv$Dm,F_surv$Dn)
  F_hull <- c(F_hull, F_hull[1])
  Equilibrium_hull <- chull(Equilibrium$Dm,Equilibrium$Dn)
  Equilibrium_hull <- c(Equilibrium_hull, Equilibrium_hull[1])
  One_hull <- chull(One_surv$Dm,One_surv$Dn)
  One_hull <- c(One_hull, One_hull[1])
  
  dat <- dat[,c(1,3,10)] ## 9 for Time to surpass, 10 for Time to Extinction and 11 for time to peak f
  
  plot(dat[,2],dat[,1],col=colors[dat[,3]], pch = 15, cex = 2.5, bty = "n", xlab = "", ylab = "")
  title(ylab = nameDm, line = 2, cex.lab = 1.5)
  mtext(nameDn, side = 1, line = 1.5, cex = 1, adj = 0.55)
  grid(lty = 2, lwd = 0.3, col = "white")
  ## Add dashed polygons
  polygon(x = HG_surv[HG_hull,2],
          y = HG_surv[HG_hull,1],
          col = adjustcolor("cyan3", alpha.f = 0.4), border = NA)
  polygon(x = F_surv[F_hull,2],
          y = F_surv[F_hull,1],
          col = adjustcolor("magenta3", alpha.f = 0.4), border = NA)
  polygon(x = Equilibrium[Equilibrium_hull,2],
          y = Equilibrium[Equilibrium_hull,1],
          col = adjustcolor("gold3", alpha.f = 0.4), border = NA)
  polygon(x = One_surv[One_hull,2],
          y = One_surv[One_hull,1],
          col = adjustcolor("lightblue3", alpha.f = 0.4), border = NA)
  
  # Overlay time to event
  points(dat[which(dat[,3]<1500),2], dat[which(dat[,3]<1500),1], col = colors[dat[which(dat[,3]<1500),3]], pch = 15, cex = 2.5)
  ## Add polygon squares
  lw <- 1.5
  polygon(x = HG_surv[HG_hull,2],
          y = HG_surv[HG_hull,1],
          col = NA, border = "cyan3", lwd = lw)
  polygon(x = F_surv[F_hull,2],
          y = F_surv[F_hull,1],
          col = NA, border = "magenta3", lwd = lw)
  polygon(x = Equilibrium[Equilibrium_hull,2],
          y = Equilibrium[Equilibrium_hull,1],
          col = NA, border = "gold3", lwd = lw)
  polygon(x = One_surv[One_hull,2],
          y = One_surv[One_hull,1],
          col = NA, border = "lightblue3", lwd = lw)
  contour(dat[,2][1:100],dat[,1][seq(1,length(dat[,1]),100)],matrix(c(dat[,3]),nrow = 100, ncol = 100), 
          col = "black", add = TRUE)
  #if (i == 3){
  #  legend(x = 0.036, y = 0.038, legend = rep("",8), pch = rep(15,8), 
  #         col = colors[seq(1,time,length.out=8)], bty = "n",  cex = 7, y.intersp = 0.15)
  #  legend(x = 0.0435, y = 0.0335, legend = leg_vals, bty = "n",  cex = 1, y.intersp = 1.25)
  #  legend(x = 0.039, y = 0.05, legend = c("Time to","extinction"), cex = 1.2, bty = "n", y.intersp = 1,
  #         text.font = 2)
  #}
  title(paste0("assimilation", " = ", e[i], " / ratio = ", ratio[i]))
}

dev.off()


################################################################################
#####################           END FIG. S3 4           ########################
################################################################################


################################################################################
#####################              FIG. S3 1            ########################
################################################################################

## Time to peak

#png("./Figures/Fig_S3_1.png", width = 1200, height = 850)
png("./Figures/Fig_S3_1.png", width = 41, height = 30, unit = "cm", res = 300)
#tiff("./Figures/Fig_S3_1.tiff", width = 1200, height = 850)

par(mfrow = c(3,3), mar = c(5,4,4,8), mai = c(0.35,0.5,0.3,0), xpd = TRUE)
for (i in 1:9){
  
  ## Load data
  dat <- readRDS(paste0("./Results/Sensitivity/Growth_r_",ratio[i],"_e_",e[i],".rds"))
  
  ## Prepare for analytical conditions
  Dm_dat <- dat$Dm
  Gm_dat <- dat$Gm
  Gn_dat <- dat$Gn
  Dn_dat <- dat$Dn
  e_dat <- dat$e
  Km_dat <- dat$Km
  Kn_dat <- dat$Kn
  
  # Lambda values for conditions
  Li <- Gm_dat-(Dm_dat*Kn_dat)
  Lii <- Gn_dat-((Dn_dat-(e_dat*Dm_dat))*Km_dat)
  
  ## Establish parameters ranges for plots
  HG_surv <- dat[which(Li > 0 & Lii < 0),c(2,4)]
  F_surv <- dat[which(Li < 0 & Lii > 0),c(2,4)]
  Equilibrium <- dat[which(Li > 0 & Lii > 0),c(2,4)]
  One_surv <- dat[which(Li < 0 & Lii < 0),c(2,4)]
  
  ## Create convex hulls with outer points for outline polygon
  HG_hull <- chull(HG_surv$Gm,HG_surv$Gn)
  HG_hull <- c(HG_hull, HG_hull[1])
  F_hull <- chull(F_surv$Gm,F_surv$Gn)
  F_hull <- c(F_hull, F_hull[1])
  Equilibrium_hull <- chull(Equilibrium$Gm,Equilibrium$Gn)
  Equilibrium_hull <- c(Equilibrium_hull, Equilibrium_hull[1])
  One_hull <- chull(One_surv$Gm,One_surv$Gn)
  One_hull <- c(One_hull, One_hull[1])
  
  dat <- dat[,c(2,4,11)] ## 9 for Time to surpass, 10 for Time to Extinction and 11 for time to peak f
  
  plot(dat[,2],dat[,1],col=colors[dat[,3]], pch = 15, cex = 2.5, bty = "n", xlab = "", ylab = "")
  title(ylab = nameGm, line = 2, cex.lab = 1.5)
  mtext(nameGn, side = 1, line = 1.5, cex = 1, adj = 0.55)
  grid(lty = 2, lwd = 0.3, col = "white")
  ## Add dashed polygons
  polygon(x = HG_surv[HG_hull,2],
          y = HG_surv[HG_hull,1],
          col = adjustcolor("cyan3", alpha.f = 0.4), border = NA)
  polygon(x = F_surv[F_hull,2],
          y = F_surv[F_hull,1],
          col = adjustcolor("magenta3", alpha.f = 0.4), border = NA)
  polygon(x = Equilibrium[Equilibrium_hull,2],
          y = Equilibrium[Equilibrium_hull,1],
          col = adjustcolor("gold3", alpha.f = 0.4), border = NA)
  polygon(x = One_surv[One_hull,2],
          y = One_surv[One_hull,1],
          col = adjustcolor("lightblue3", alpha.f = 0.4), border = NA)
  
  # Overlay time to event
  points(dat[which(dat[,3]<1500),2], dat[which(dat[,3]<1500),1], col = colors[dat[which(dat[,3]<1500),3]], pch = 15, cex = 2.5)
  ## Add polygon squares
  lw <- 1.5
  polygon(x = HG_surv[HG_hull,2],
          y = HG_surv[HG_hull,1],
          col = NA, border = "cyan3", lwd = lw)
  polygon(x = F_surv[F_hull,2],
          y = F_surv[F_hull,1],
          col = NA, border = "magenta3", lwd = lw)
  polygon(x = Equilibrium[Equilibrium_hull,2],
          y = Equilibrium[Equilibrium_hull,1],
          col = NA, border = "gold3", lwd = lw)
  polygon(x = One_surv[One_hull,2],
          y = One_surv[One_hull,1],
          col = NA, border = "lightblue3", lwd = lw)
  contour(dat[,2][1:100],dat[,1][seq(1,length(dat[,1]),100)],matrix(c(dat[,3]),nrow = 100, ncol = 100), 
          col = "black", add = TRUE)
  i#f (i == 3){
  #  legend(x = 0.0645, y = 0.019, legend = rep("",8), pch = rep(15,8), 
  #         col = colors[seq(1,time,length.out=8)], bty = "n",  cex = 7, y.intersp = 0.15)
  #  legend(x = 0.075, y = 0.0165, legend = leg_vals, bty = "n",  cex = 1, y.intersp = 1.25)
  #  legend(x = 0.0685, y = 0.02, legend = c("Time to","peak"), cex = 1.2, bty = "n", y.intersp = 1,
  #         text.font = 2)
  #}
  title(paste0("assimilation", " = ", e[i], " / ratio = ", ratio[i]))
}

dev.off()


################################################################################
#####################           END FIG. S3 1           ########################
################################################################################


################################################################################
#####################              FIG. S3 2            ########################
################################################################################

## Time to surpass

#png("./Figures/Fig_S3_2.png", width = 1200, height = 850)
png("./Figures/Fig_S3_2.png", width = 41, height = 30, unit = "cm", res = 300)
#tiff("./Figures/Fig_S3_2.tiff", width = 1200, height = 850)

par(mfrow = c(3,3), mar = c(5,4,4,8), mai = c(0.35,0.5,0.3,0), xpd = TRUE)
for (i in 1:9){
  
  ## Load data
  dat <- readRDS(paste0("./Results/Sensitivity/Growth_r_",ratio[i],"_e_",e[i],".rds"))
  
  ## Prepare for analytical conditions
  Dm_dat <- dat$Dm
  Gm_dat <- dat$Gm
  Gn_dat <- dat$Gn
  Dn_dat <- dat$Dn
  e_dat <- dat$e
  Km_dat <- dat$Km
  Kn_dat <- dat$Kn
  
  # Lambda values for conditions
  Li <- Gm_dat-(Dm_dat*Kn_dat)
  Lii <- Gn_dat-((Dn_dat-(e_dat*Dm_dat))*Km_dat)
  
  ## Establish parameters ranges for plots
  HG_surv <- dat[which(Li > 0 & Lii < 0),c(2,4)]
  F_surv <- dat[which(Li < 0 & Lii > 0),c(2,4)]
  Equilibrium <- dat[which(Li > 0 & Lii > 0),c(2,4)]
  One_surv <- dat[which(Li < 0 & Lii < 0),c(2,4)]
  
  ## Create convex hulls with outer points for outline polygon
  HG_hull <- chull(HG_surv$Gm,HG_surv$Gn)
  HG_hull <- c(HG_hull, HG_hull[1])
  F_hull <- chull(F_surv$Gm,F_surv$Gn)
  F_hull <- c(F_hull, F_hull[1])
  Equilibrium_hull <- chull(Equilibrium$Gm,Equilibrium$Gn)
  Equilibrium_hull <- c(Equilibrium_hull, Equilibrium_hull[1])
  One_hull <- chull(One_surv$Gm,One_surv$Gn)
  One_hull <- c(One_hull, One_hull[1])
  
  dat <- dat[,c(2,4,9)] ## 9 for Time to surpass, 10 for Time to Extinction and 11 for time to peak f
  
  plot(dat[,2],dat[,1],col=colors[dat[,3]], pch = 15, cex = 2.5, bty = "n", xlab = "", ylab = "")
  title(ylab = nameGm, line = 2, cex.lab = 1.5)
  mtext(nameGn, side = 1, line = 1.5, cex = 1, adj = 0.55)
  grid(lty = 2, lwd = 0.3, col = "white")
  ## Add dashed polygons
  polygon(x = HG_surv[HG_hull,2],
          y = HG_surv[HG_hull,1],
          col = adjustcolor("cyan3", alpha.f = 0.4), border = NA)
  polygon(x = F_surv[F_hull,2],
          y = F_surv[F_hull,1],
          col = adjustcolor("magenta3", alpha.f = 0.4), border = NA)
  polygon(x = Equilibrium[Equilibrium_hull,2],
          y = Equilibrium[Equilibrium_hull,1],
          col = adjustcolor("gold3", alpha.f = 0.4), border = NA)
  polygon(x = One_surv[One_hull,2],
          y = One_surv[One_hull,1],
          col = adjustcolor("lightblue3", alpha.f = 0.4), border = NA)
  
  # Overlay time to event
  points(dat[which(dat[,3]<1500),2], dat[which(dat[,3]<1500),1], col = colors[dat[which(dat[,3]<1500),3]], pch = 15, cex = 2.5)
  ## Add polygon squares
  lw <- 1.5
  polygon(x = HG_surv[HG_hull,2],
          y = HG_surv[HG_hull,1],
          col = NA, border = "cyan3", lwd = lw)
  polygon(x = F_surv[F_hull,2],
          y = F_surv[F_hull,1],
          col = NA, border = "magenta3", lwd = lw)
  polygon(x = Equilibrium[Equilibrium_hull,2],
          y = Equilibrium[Equilibrium_hull,1],
          col = NA, border = "gold3", lwd = lw)
  polygon(x = One_surv[One_hull,2],
          y = One_surv[One_hull,1],
          col = NA, border = "lightblue3", lwd = lw)
  contour(dat[,2][1:100],dat[,1][seq(1,length(dat[,1]),100)],matrix(c(dat[,3]),nrow = 100, ncol = 100), 
          col = "black", add = TRUE)
  #if (i == 3){legend(x = 0.0645, y = 0.019, legend = rep("",8), pch = rep(15,8), 
  #                   col = colors[seq(1,time,length.out=8)], bty = "n",  cex = 7, y.intersp = 0.15)
  #  legend(x = 0.075, y = 0.0165, legend = leg_vals, bty = "n",  cex = 1, y.intersp = 1.25)
  #  legend(x = 0.0685, y = 0.02, legend = c("Time to","surpass"), cex = 1.2, bty = "n", y.intersp = 1,
  #         text.font = 2)}
  
  title(paste0("assimilation", " = ", e[i], " / ratio = ", ratio[i]))
}

dev.off()

################################################################################
#####################           END FIG. S3 2           ########################
################################################################################


