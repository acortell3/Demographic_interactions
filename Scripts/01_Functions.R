
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
################                   SECTION 1                    ################
################                                                ################
################################################################################

## This section contains the different packages and functions needed for later sections.

################################################################################
################           LOAD PACKAGES AND FUNCTIONS          ################
################################################################################


## Load libraries
library(rcarbon)  ## For radiocarbon calibration and SPDs
library(stats)    ## For distributions
library(foreach)  ## For parallel computing
library(doParallel) ## For parallel computing
library(deSolve) ## For differential equations

## For plots
library(coda) ## HDP
library(colorRamps) ## Colors


## FUNCTIONS
## Function 1. LotVol function
LotVolme <- function(pars, times=seq(0, 1000, by = 1)) {
  derivs <- function(t, state, pars) { # returns rate of change
    with(as.list(c(state, pars)), {
      Im <- Dm * M * N ## Interaction effect
      H <- Gm * M * (1 - M/Km) ## Hg pop increase
      Phi <- (Gn + m) * N * (1 - N/Kn) ## F pop increase
      In <- Dn * M * N
      
      dM    <- H - Im
      dN    <- Im * e + Phi - In
      
      return(list(c(dM, dN)))
    })
  }
  
  return(ode(y = state, times = times, func = derivs, parms = pars))
}

## Function 2. Extracts dates according to the probability that one site exists and it is meso/neo
cult_afil <- function(x,y,start,end,dates){
  
  ## x is meso population
  ## y is neo population
  ## start is year where the model starts
  ## end is year where the model ends
  ## dates is how many dates do we want
  
  x[x < 0] <- 0
  o <- x+y ## Total population
  pm <- x/o ## Meso probability
  
  
  ts <- seq(start,end) ## Time span
  
  ye <- sample(ts,dates,prob=o) ## Number of dates selected
  
  m <- rbinom(length(ye),1,pm[ts%in%ye]) ## Are the dates Mesolithic?
  
  m[is.na(m)] <- 1
  
  ## Solution for the case where there are no mesolithic dates
  if (sum(m)<1){
    m[1] <- 1
  }
  
  if (sum(m) == length(m)){
    m[length(m)] <- 0
  }
  
  r <- data.frame("Year" = sort(ye, decreasing = TRUE),
                  "Culture" = as.factor(ifelse(m==1,"hg","f"))) ## output
  
  return(r)
}

#res <- cult_afil(x = dat[,1], y = dat[,2], dates = 100, start = 6600, end = 5600)

## Examples

#par(mfrow = c(1,2))
#plot(dat[,1], type = "l", ylim = c(0,100), col = "darkred", lwd = 1.5, xlab = "Pop", ylab = "Time")
#lines(dat[,2], col = "darkblue", lwd = 1.5)

#plot(res, xlim = rev(range(res$Year)), pch = 16, col = "blue", yaxt = "n")
#axis(2, at=1:2, labels = c("f","hg"), font = 2, las = 1)
#res

## Function 3. Select dates whose probability density is above/below a specific threshold
## Requires Rcarbon
## x = the date [0,inf]
## s = the standard deviation [0,inf]
## ot = older threshold [0,inf]
## rt = recent threshold [0,inf]
## return = How much density do we want in the thershold [0,1]

Date_pres <- function(x,s,ot,rt){
  cd <- calibrate(x,s)
  pr <- cd$grids$`1`[cd$grids$`1`$calBP<ot&cd$grids$`1`$calBP>rt,]
  return(sum(pr$PrDens))
}


##### FUNCTIONS FOR PLOTS
## Functions for computing max and min values
maxmin <- function(x){
  maxval <- max(x)
  minval <- min(x)
  return(c(minval,maxval))
}
