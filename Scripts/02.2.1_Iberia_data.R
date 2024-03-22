
################################################################################
################################################################################

####### Code for: Demographic interactions of the last hunter-gatherers and the first farmers (2024):


#### License: Permission is granted to use and adapt this code. Please 
#### acknowledge authorship when appropriate.

################################################################################
################                                                ################
################                  SECTION 2.2.1                 ################
################                                                ################
################################################################################

### OBSERVED IBERIA

## Load data and libraries
library(readxl)
library(rcarbon)

set.seed(123)

### DATES

Dat_dat <- read.csv("Data/Observed/Iberia_selected_ogp.csv")

All <- Dat_dat[,c(1,6,3,4,5)] ## Need to add region
All$Fecha.BP <- as.numeric(All$Fecha.BP)
All$Desviacion <- as.numeric(All$Desviacion)
All <- na.omit(All)
All <- All[All$Desviacion <= 120,]

## EXPLORATION OF THE DATASET

## Iberia
Neo_Iberia_dates <- All[All$Adscripcion == "Neolithic",]
Meso_Iberia_dates <- All[All$Adscripcion == "Mesolithic",]

## We prepare the bins. We will use 50 years
Iberia_hg_bin <- binPrep(sites=Meso_Iberia_dates$Yacimiento, ages = Meso_Iberia_dates$Fecha.BP, h = 20)
Iberia_f_bin <- binPrep(sites=Neo_Iberia_dates$Yacimiento, ages = Neo_Iberia_dates$Fecha.BP, h = 20)

## Use random thinning on the bins. Maximum of five dates per site
## Calibrate dates
Meso_cal_Iberia <- calibrate(Meso_Iberia_dates$Fecha.BP, Meso_Iberia_dates$Desviacion,calCurves = 'intcal20')
Neo_cal_Iberia <- calibrate(Neo_Iberia_dates$Fecha.BP, Neo_Iberia_dates$Desviacion,calCurves = 'intcal20')

Iberia_hg_thin <- Meso_cal_Iberia[thinDates(ages = Meso_Iberia_dates$Fecha.BP, errors = Meso_Iberia_dates$Desviacion, 
                                          bins = Iberia_hg_bin, size = 1, method = "random")]

Iberia_f_thin <- Neo_cal_Iberia[thinDates(ages = Neo_Iberia_dates$Fecha.BP, errors = Neo_Iberia_dates$Desviacion, 
                                        bins = Iberia_f_bin, size = 1, method = "random")]

Neo_spd_Iberia <- spd(Iberia_f_thin, timeRange = c(8200,7300))
Meso_spd_Iberia <- spd(Iberia_hg_thin, timeRange = c(8200,7300))

Total_dates_Iberia <- Neo_spd_Iberia$metadata$ndates + Meso_spd_Iberia$metadata$ndates
Total_dates_Iberia

par(mfrow = c(1,1))
plot(Neo_spd_Iberia, main = paste0("Iberia dates (ndates = ",Total_dates_Iberia," )"))
plot(Meso_spd_Iberia, add = TRUE)
plot(Neo_spd_Iberia, type = "simple", col = "red", add = TRUE)
plot(Meso_spd_Iberia, type = "simple", col = "blue", add = TRUE)

### EXPORT DATASET

## To store different beginning
beginning <- seq(100,400,100)

Newest <- 7300 ## Based on observation of the SPD
Oldest <- c(Newest + 352 + beginning[1], ## 750 is based on the observation of Neo_spd_Iberia$grid$PrDens
            Newest + 352 + beginning[2],
            Newest + 352 + beginning[3],
            Newest + 352 + beginning[4],
            Newest + 352 + beginning[5])

## Generate dates with different beginnings

for (j in 1:length(beginning)){
  ## Select dates only if an amount of their probability is within the time period of study
  oldest_Dat <- Oldest[j]
  recent_Dat <- Newest
  
  Dat_hg_dates_sel <- Meso_Iberia_dates
  
  for (i in 1:nrow(Dat_hg_dates_sel)){
    cond <- Date_pres(x=Meso_Iberia_dates[i,3],s=Meso_Iberia_dates[i,4],ot=oldest_Dat,rt=recent_Dat)
    if (cond >= 0.5){ ## 50% of the probability accepted
      Dat_hg_dates_sel[i,3] <- Meso_Iberia_dates[i,3]
      Dat_hg_dates_sel[i,4] <- Meso_Iberia_dates[i,4]
    } else {
      Dat_hg_dates_sel[i,3] <- NA
      Dat_hg_dates_sel[i,4] <- NA
    }
  }
  
  Dat_hg_dates_sel <- Dat_hg_dates_sel[complete.cases(Dat_hg_dates_sel),]
  
  Dat_f_dates_sel <- Neo_Iberia_dates
  
  for (i in 1:nrow(Dat_f_dates_sel)){
    cond <- Date_pres(x=Neo_Iberia_dates[i,3],s=Neo_Iberia_dates[i,4],ot=oldest_Dat,rt=recent_Dat)
    if (cond >= 0.5){ ## 50% of the probability accepted
      Dat_f_dates_sel[i,3] <- Neo_Iberia_dates[i,3]
      Dat_f_dates_sel[i,4] <- Neo_Iberia_dates[i,4]
    } else {
      Dat_f_dates_sel[i,3] <- NA
      Dat_f_dates_sel[i,4] <- NA
    }
  }
  
  Dat_f_dates_sel <- Dat_f_dates_sel[complete.cases(Dat_f_dates_sel),]
  
  
  # Produce the spds with the selected dates
  ## bins
  Dat_hg_bin_sel <- binPrep(sites=Dat_hg_dates_sel$Yacimiento, ages = Dat_hg_dates_sel$Fecha.BP, h = 20)
  Dat_f_bin_sel <- binPrep(sites=Dat_f_dates_sel$Yacimiento, ages = Dat_f_dates_sel$Fecha.BP, h = 20)
  
  Dat_hg_cal_dates_sel <- calibrate(Dat_hg_dates_sel$Fecha.BP,Dat_hg_dates_sel$Desviacion,calCurves = 'intcal20')
  Dat_f_cal_dates_sel <- calibrate(Dat_f_dates_sel$Fecha.BP,Dat_f_dates_sel$Desviacion,calCurves = 'intcal20')
  
  Dat_hg_thin_sel <- Dat_hg_cal_dates_sel[thinDates(ages = Dat_hg_dates_sel$Fecha.BP, errors = Dat_hg_dates_sel$Desviacion, 
                                                    bins = Dat_hg_bin_sel, size = 1, method = "random")]
  
  Dat_f_thin_sel <- Dat_f_cal_dates_sel[thinDates(ages = Dat_f_dates_sel$Fecha.BP, errors = Dat_f_dates_sel$Desviacion, 
                                                  bins = Dat_f_bin_sel, size = 1, method = "random")]
  
  Dat_hg_spd_sel <- spd(Dat_hg_thin_sel, timeRange = c(oldest_Dat,recent_Dat))
  Dat_f_spd_sel <- spd(Dat_f_thin_sel, timeRange = c(oldest_Dat,recent_Dat))
  
  Total_dates_Iberia <- Dat_f_spd_sel$metadata$ndates + Dat_hg_spd_sel$metadata$ndates
  Total_dates_Iberia
  
  par(mfrow = c(1,1))
  plot(Dat_f_spd_sel, main = paste0("Iberia dates (ndates = ",Total_dates_Iberia," )"))
  plot(Dat_hg_spd_sel, add = TRUE)
  plot(Dat_f_spd_sel, type = "simple", col = "red", add = TRUE)
  plot(Dat_hg_spd_sel, type = "simple", col = "blue", add = TRUE)
  
  saveRDS(Dat_hg_spd_sel,file = paste0("Data/Observed/Iberia/Iberia_",as.character(beginning[j]),"_years_before_hg.rds"))
  saveRDS(Dat_f_spd_sel,file = paste0("Data/Observed/Iberia/Iberia_",as.character(beginning[j]),"_years_before_f.rds"))
  saveRDS(Dat_hg_dates_sel$Desviacion, file = paste0("Data/Observed/Iberia/Iberia_",as.character(beginning[j]),"_years_before_Sd_hg.rds"))
  saveRDS(Dat_f_dates_sel$Desviacion, file = paste0("Data/Observed/Iberia/Iberia_",as.character(beginning[j]),"_years_before_Sd_f.rds"))
}



