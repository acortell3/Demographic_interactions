
### OBSERVED DENMARK

## Load data and libraries
library(readxl)
library(rcarbon)

set.seed(123)

### DATES

Dat_dat <- read.csv2("Data/Observed/Scandinavia_Meso_Neo_dates2.csv")

All <- Dat_dat[,c(2,3,6,7,4)]
All$C14.Age <- as.numeric(All$C14.Age)
All$C14.SD <- as.numeric(All$C14.SD)
All <- na.omit(All)
All <- All[All$C14.SD <= 120,]

## EXPLORATION OF THE DATASET

## Denmark
Neo_Denmark_dates <- All[All$Culture1 == "Trichterbecher" & All$Country == "Denmark",]
Meso_Denmark_dates <- All[All$Culture1 == "Ertebølle" & All$Country == "Denmark",]

## Try removing potential outliers
Neo_Denmark_dates <- Neo_Denmark_dates[Neo_Denmark_dates$C14.Age < 5700,]

## We prepare the bins. We will use 50 years
Denmark_hg_bin <- binPrep(sites=Meso_Denmark_dates$SiteID, ages = Meso_Denmark_dates$C14.Age, h = 20)
Denmark_f_bin <- binPrep(sites=Neo_Denmark_dates$SiteID, ages = Neo_Denmark_dates$C14.Age, h = 20)

## Use random thinning on the bins. Maximum of five dates per site
## Calibrate dates
Meso_cal_Denmark <- calibrate(Meso_Denmark_dates$C14.Age, Meso_Denmark_dates$C14.SD,calCurves = 'intcal20')
Neo_cal_Denmark <- calibrate(Neo_Denmark_dates$C14.Age, Neo_Denmark_dates$C14.SD,calCurves = 'intcal20')

Denmark_hg_thin <- Meso_cal_Denmark[thinDates(ages = Meso_Denmark_dates$C14.Age, errors = Meso_Denmark_dates$C14.SD, 
                                          bins = Denmark_hg_bin, size = 1, method = "random")]

Denmark_f_thin <- Neo_cal_Denmark[thinDates(ages = Neo_Denmark_dates$C14.Age, errors = Neo_Denmark_dates$C14.SD, 
                                        bins = Denmark_f_bin, size = 1, method = "random")]

Neo_spd_Denmark <- spd(Denmark_f_thin, timeRange = c(7000,5000))
Meso_spd_Denmark <- spd(Denmark_hg_thin, timeRange = c(7000,5000))

Total_dates_Denmark <- Neo_spd_Denmark$metadata$ndates + Meso_spd_Denmark$metadata$ndates
Total_dates_Denmark

par(mfrow = c(1,1))
plot(Neo_spd_Denmark, main = paste0("Denmark dates (ndates = ",Total_dates_Denmark," )"))
plot(Meso_spd_Denmark, add = TRUE)
plot(Neo_spd_Denmark, type = "simple", col = "red", add = TRUE)
plot(Meso_spd_Denmark, type = "simple", col = "blue", add = TRUE)

### EXPORT DATASET

## To store different beginning
beginning <- seq(100,400,100)

Newest <- 5600 ## Based on observation of the SPD
Oldest <- c(Newest + 695 + beginning[1], ## 695 is based on the observation of the SPD
            Newest + 695 + beginning[2],
            Newest + 695 + beginning[3],
            Newest + 695 + beginning[4],
            Newest + 695 + beginning[5])

## Generate dates with different beginnings

for (j in 1:length(beginning)){
  ## Select dates only if an amount of their probability is within the time period of study
  oldest_Dat <- Oldest[j]
  recent_Dat <- Newest
  
  Dat_hg_dates_sel <- Meso_Denmark_dates
  
  for (i in 1:nrow(Dat_hg_dates_sel)){
    cond <- Date_pres(x=Meso_Denmark_dates[i,3],s=Meso_Denmark_dates[i,4],ot=oldest_Dat,rt=recent_Dat)
    if (cond >= 0.5){ ## 50% of the probability accepted
      Dat_hg_dates_sel[i,3] <- Meso_Denmark_dates[i,3]
      Dat_hg_dates_sel[i,4] <- Meso_Denmark_dates[i,4]
    } else {
      Dat_hg_dates_sel[i,3] <- NA
      Dat_hg_dates_sel[i,4] <- NA
    }
  }
  
  Dat_hg_dates_sel <- Dat_hg_dates_sel[complete.cases(Dat_hg_dates_sel),]
  
  Dat_f_dates_sel <- Neo_Denmark_dates
  
  for (i in 1:nrow(Dat_f_dates_sel)){
    cond <- Date_pres(x=Neo_Denmark_dates[i,3],s=Neo_Denmark_dates[i,4],ot=oldest_Dat,rt=recent_Dat)
    if (cond >= 0.5){ ## 50% of the probability accepted
      Dat_f_dates_sel[i,3] <- Neo_Denmark_dates[i,3]
      Dat_f_dates_sel[i,4] <- Neo_Denmark_dates[i,4]
    } else {
      Dat_f_dates_sel[i,3] <- NA
      Dat_f_dates_sel[i,4] <- NA
    }
  }
  
  Dat_f_dates_sel <- Dat_f_dates_sel[complete.cases(Dat_f_dates_sel),]
  
  
  # Produce the spds with the selected dates
  ## bins
  Dat_hg_bin_sel <- binPrep(sites=Dat_hg_dates_sel$SiteID, ages = Dat_hg_dates_sel$C14.Age, h = 20)
  Dat_f_bin_sel <- binPrep(sites=Dat_f_dates_sel$SiteID, ages = Dat_f_dates_sel$C14.Age, h = 20)
  
  Dat_hg_cal_dates_sel <- calibrate(Dat_hg_dates_sel$C14.Age,Dat_hg_dates_sel$C14.SD,calCurves = 'intcal20')
  Dat_f_cal_dates_sel <- calibrate(Dat_f_dates_sel$C14.Age,Dat_f_dates_sel$C14.SD,calCurves = 'intcal20')
  
  Dat_hg_thin_sel <- Dat_hg_cal_dates_sel[thinDates(ages = Dat_hg_dates_sel$C14.Age, errors = Dat_hg_dates_sel$C14.SD, 
                                                    bins = Dat_hg_bin_sel, size = 1, method = "random")]
  
  Dat_f_thin_sel <- Dat_f_cal_dates_sel[thinDates(ages = Dat_f_dates_sel$C14.Age, errors = Dat_f_dates_sel$C14.SD, 
                                                  bins = Dat_f_bin_sel, size = 1, method = "random")]
  
  Dat_hg_spd_sel <- spd(Dat_hg_thin_sel, timeRange = c(oldest_Dat,recent_Dat))
  Dat_f_spd_sel <- spd(Dat_f_thin_sel, timeRange = c(oldest_Dat,recent_Dat))
  
  Total_dates_Denmark <- Dat_f_spd_sel$metadata$ndates + Dat_hg_spd_sel$metadata$ndates
  Total_dates_Denmark
  
  par(mfrow = c(1,1))
  plot(Dat_f_spd_sel, main = paste0("Denmark dates (ndates = ",Total_dates_Denmark," )"))
  plot(Dat_hg_spd_sel, add = TRUE)
  plot(Dat_f_spd_sel, type = "simple", col = "red", add = TRUE)
  plot(Dat_hg_spd_sel, type = "simple", col = "blue", add = TRUE)
  
  saveRDS(Dat_hg_spd_sel,file = paste0("Data/Observed/Denmark/Denmark_",as.character(beginning[j]),"_years_before_hg.rds"))
  saveRDS(Dat_f_spd_sel,file = paste0("Data/Observed/Denmark/Denmark_",as.character(beginning[j]),"_years_before_f.rds"))
  saveRDS(Dat_hg_dates_sel$C14.SD, file = paste0("Data/Observed/Denmark/Denmark_",as.character(beginning[j]),"_years_before_Sd_hg.rds"))
  saveRDS(Dat_f_dates_sel$C14.SD, file = paste0("Data/Observed/Denmark/Denmark_",as.character(beginning[j]),"_years_before_Sd_f.rds"))
}



