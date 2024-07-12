
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
################                   SECTION 54                   ################
################                                                ################
################################################################################

################################################################################
########################     SENSITIVITY ANALYSIS     ##########################
################################################################################

source("Scripts/01_Functions.R")

## Assessment of growth parameters

## Utilities
le <- 100 ## length of explored parameter
time <- 1500 ## length of the process

################################################################################
########################    GROWTH(s) PARAMETER(s)    ##########################
################################################################################

pars_used <- data.frame("Dm" = rep(NA,le),
                        "Gm" = rep(NA,le),
                        "Dn" = rep(NA,le),
                        "Gn" = rep(NA,le),
                        "e" = rep(NA,le),
                        "Km" = rep(NA,le),
                        "Kn" = rep(NA,le),
                        "initial.ratio" = rep(NA,le),
                        "time.to.overcome" = rep(NA,le),
                        "time.to.extinction.hg" = rep(NA,le),
                        "time.to.peak.f" = rep(NA,le))

final_pars <- data.frame()
pop_growth <- list()

## Parameters
e <- c(0,0.1,0.2) ## Interaction
Dm <- 0.01 ## Mortality hg
Dn <- 0.01 ## Mortality f
Gm <- seq(from = 0.001, to = 0.022, length.out = le)
Gn <- seq(from = 0.015, to = 0.07, length.out = le)

M <- 1 ## Initial population hg
N <- c(0.1,0.2,0.4) ## Initial population f

for (h in 1:3){ ## Set population ratio
  state <- c(M = M, N = N[h])
  
  for (k in 1:3){ ## Set interaction parameter
    
    for (i in 1:le){
      ## It runs all the Gn on a single Gm and then it goes to the next Gm and does the same
      Gm_c <- Gm[i]
      for (j in 1:le){
        ## Prepare data frame
        Dm_c <- Dm
        Dn_c <- Dn
        Gn_c <- Gn[j]
        e_c <- e[k]
        Km_c <- M
        Kn_c <- M * 1.5
        
        
        pars <- list("Dm" = Dm_c,
                     "Gm" = Gm_c,
                     "Dn" = Dn_c,
                     "Gn" = Gn_c,
                     "e" = e_c,
                     "Km" = Km_c,
                     "Kn" = Kn_c,
                     "m" = 0) ## m is assumed in Gn
        
        
        sim <- LotVolme(pars = pars, times = c(1:time))
        
        ratio <- N[h]
        time.to.overcome <- which(sim[,2] > sim[,3])[length(which(sim[,2] > sim [,3]))]+1
        time.to.overcome[is.na(time.to.overcome)] <- time ## If there is no overcome
        time.to.extinction.hg <- which(sim[,2] < 0.01)[1] ## It will never be 0, so we choose 0.001
        time.to.extinction.hg[is.na(time.to.extinction.hg)] <- time ## If there is no extinction
        time.to.peak.f <- which(sim[,3] >= Kn_c - 0.01)[1]
        time.to.peak.f[is.na(time.to.peak.f)] <- time ## If there is no farmer peak after overcome
        
        vec_to_pass <- c(Dm_c,Gm_c,Dn_c,Gn_c,e_c,Km_c,Kn_c,ratio,time.to.overcome,time.to.extinction.hg,time.to.peak.f)
        
        pop_growth[[length(pop_growth) + 1]] <- sim
        pars_used[j,] <- vec_to_pass
      }
      final_pars <- rbind(final_pars,pars_used)
      
    }
    
    saveRDS(final_pars,paste0("Results/Sensitivity/Growth_r_",ratio,"_e_",e_c,".rds"))
    saveRDS(pop_growth,paste0("Results/Sensitivity/Growth_pop_r_",ratio,"_e_",e_c,".rds"))
    final_pars <- data.frame()
    pop_growth <- list()
    
  }
    
}
  

################################################################################
########################     DEATH(s) PARAMETER(s)    ##########################
################################################################################

pars_used <- data.frame("Dm" = rep(NA,le),
                        "Gm" = rep(NA,le),
                        "Dn" = rep(NA,le),
                        "Gn" = rep(NA,le),
                        "e" = rep(NA,le),
                        "Km" = rep(NA,le),
                        "Kn" = rep(NA,le),
                        "initial.ratio" = rep(NA,le),
                        "time.to.overcome" = rep(NA,le),
                        "time.to.extinction.hg" = rep(NA,le),
                        "time.to.peak.f" = rep(NA,le))

final_pars <- data.frame()
pop_death <- list()

## Parameters
e <- c(0,0.1,0.2) ## Interaction
Dm <- seq(from = 0.001, to = 0.04, length.out = le)
Dn <- seq(from = 0.001, to = 0.04, length.out = le)
Gm <- 0.015
Gn <- 0.02 ## Just a little bit more farmers

M <- 1 ## Initial population hg
N <- c(0.1,0.2,0.4) ## Initial population f

for (h in 1:3){ ## Set population ratio
  state <- c(M = M, N = N[h])
  
  for (k in 1:3){ ## Set interaction parameter
    
    for (i in 1:le){
      ## It runs all the Gn on a single Gm and then it goes to the next Gm and does the same
      Dm_c <- Dm[i]
      for (j in 1:le){
        ## Prepare data frame
        Gm_c <- Gm
        Gn_c <- Gn
        Dn_c <- Dn[j]
        e_c <- e[k]
        Km_c <- M
        Kn_c <- M * 1.5
        
        
        pars <- list("Dm" = Dm_c,
                     "Gm" = Gm_c,
                     "Dn" = Dn_c,
                     "Gn" = Gn_c,
                     "e" = e_c,
                     "Km" = Km_c,
                     "Kn" = Kn_c,
                     "m" = 0) ## m is assumed in Gn
        
        
        sim <- LotVolme(pars = pars, times = c(1:time))
        
        ratio <- N[h]
        time.to.overcome <- which(sim[,2] > sim[,3])[length(which(sim[,2] > sim [,3]))]+1
        time.to.overcome[is.na(time.to.overcome)] <- time ## If there is no overcome
        time.to.extinction.hg <- which(sim[,2] < 0.01)[1] ## It will never be 0, so we choose 0.001
        time.to.extinction.hg[is.na(time.to.extinction.hg)] <- time ## If there is no extinction
        time.to.peak.f <- which(sim[,3] >= Kn_c - 0.01)[1]
        time.to.peak.f[is.na(time.to.peak.f)] <- time ## If there is no farmer peak after overcome
        
        vec_to_pass <- c(Dm_c,Gm_c,Dn_c,Gn_c,e_c,Km_c,Kn_c,ratio,time.to.overcome,time.to.extinction.hg,time.to.peak.f)
        
        pop_death[[length(pop_death) + 1]] <- sim
        pars_used[j,] <- vec_to_pass
      }
      final_pars <- rbind(final_pars,pars_used)
    }
    
    saveRDS(final_pars,paste0("Results/Sensitivity/Death_r_",ratio,"_e_",e_c,".rds"))
    saveRDS(pop_death,paste0("Results/Sensitivity/Death_pop_r_",ratio,"_e_",e_c,".rds"))
    final_pars <- data.frame()
    pop_death <- list()
  }
  
}

