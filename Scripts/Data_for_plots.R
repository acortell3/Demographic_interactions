

### DATA FOR PLOTS

## Load selected simulations

# Simulated
model_list_sim <- readRDS("./Data/Sim_selected_500/model_sim_best500.rds")
fit_list_sim <- readRDS("./Data/Sim_selected_500/fit_sim_best500.rds")

# Iberia
model_list_Ib <- readRDS("./Data/Observed/Iberia/Best_500/model_Iberia_best500.rds")
fit_list_Ib <- readRDS("./Data/Observed/Iberia/Best_500/fit_Iberia_best500.rds")

# Japan
model_list_Ja <- readRDS("./Data/Observed/Japan/Best_500/model_Japan_best500.rds")
fit_list_Ja <- readRDS("./Data/Observed/Japan/Best_500/fit_Japan_best500.rds")

# Denmark
model_list_De <- readRDS("./Data/Observed/Denmark/Best_500/model_Denmark_best500.rds")
fit_list_De <- readRDS("./Data/Observed/Denmark/Best_500/fit_Denmark_best500.rds")

## Load observed SPDs

# Simulated
Neo_spd_sim <- readRDS("Data/Simulated/1000_years/Tactical_100_dates_f_1000_years_2.rds")
Meso_spd_sim <- readRDS("Data/Simulated/1000_years/Tactical_100_dates_hg_1000_years_2.rds")
Sd_hg_sim <- readRDS("Data/Simulated/1000_years/Sd_Tactical_100_dates_hg_1000_years_2.rds")
Sd_f_sim <- readRDS("Data/Simulated/1000_years/Sd_Tactical_100_dates_f_1000_years_2.rds")

# Iberia
Neo_spd_Ib <- readRDS("Data/Observed/Iberia/Iberia_100_years_before_f.rds")
Meso_spd_Ib <- readRDS("Data/Observed/Iberia/Iberia_100_years_before_hg.rds")
Sd_hg_Ib <- readRDS("Data/Observed/Iberia/Iberia_100_years_before_Sd_hg.rds")
Sd_f_Ib <- readRDS("Data/Observed/Iberia/Iberia_100_years_before_Sd_f.rds")

# Japan
Neo_spd_Ja <- readRDS("Data/Observed/Japan/Japan_100_years_before_f.rds")
Meso_spd_Ja <- readRDS("Data/Observed/Japan/Japan_100_years_before_hg.rds")
Sd_hg_Ja <- readRDS("Data/Observed/Japan/Japan_100_years_before_Sd_hg.rds")
Sd_f_Ja <- readRDS("Data/Observed/Japan/Japan_100_years_before_Sd_f.rds")

# Denmark
Neo_spd_De <- readRDS("Data/Observed/Denmark/Denmark_100_years_before_f.rds")
Meso_spd_De <- readRDS("Data/Observed/Denmark/Denmark_100_years_before_hg.rds")
Sd_hg_De <- readRDS("Data/Observed/Denmark/Denmark_100_years_before_Sd_hg.rds")
Sd_f_De <- readRDS("Data/Observed/Denmark/Denmark_100_years_before_Sd_f.rds")







