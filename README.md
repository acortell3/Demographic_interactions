

<!-- rmarkdown v1 -->

Readme for the reproducibility of the paper "Demographic Interactions between the last hunter-gatherers and the first farmers", by A. Cortell-Nicolau, J. Rivas, E. R. Crema, S. Shennan, O. García-Puchol, J. Kolár, R. Staniuk and A. Timpson

## STRUCTURE

The present repository contains the following four folders:

### Data

Due to the large volume of data generated (~36Gb), here we only present the data necessary to run the scripts (e.g. the results).

This is:

- `csv` files to create the SPDs.
- SPDs both the simulated ones and the SPDs per region.
- Best 500 simulations for each dataset

Please note that for the tactical simulation, the different exploratory analyses (e.g. tries with different starts and different lengths of the process) are not included. Only the ones used and shown in the manuscript.

The folder `Data` contains three subfolders: `Sim_selected_500`, `Simulated` and `Observed`. The first one contains only the best 500 fits for the simulated data and the second one the information for the simulated SPDs. As for the folder `Observed`, it contains the `.csv` with the dates of each dataset and three more subfolders (`Denmark`, `Iberia` and `Japan`) each of which contains the SPDs and another subfolder with the best 500 fits. If the structure of this directory is not changed, all scripts should work without the need to reroute paths. 

Additionally, the full output of the results has been uploaded to Zenodo for accountabiliity (doi: 10.5281/zenodo.10837744) although, as mentioned before, it is not necessary to reproduce the code (it's actually the results that the code produces!).


### Figures

Contains the figures used in the paper. These can be generated again using the script `05_Plots_Results.R`. **In order to reproduce figure 6, it is necessary to run the script `04_Sensitivity.R` first**.


### Scripts

The scripts used for the paper. These can be conceived as follows:

`01_Functions.R`
Contains all the necessary custom functions and (most) necessary packages to run the rest of scripts

`02*_data.R`
These generate the simulated data and prepare and filter the observed data according to the criteria specified in the main manuscript.

`03*ABC_*.R`
The scripts where the model is done. They retrieve the data from the structured folders and produced in the scripts above.

`04_Sensitivity.R`
Calculations to explore the parameter space. **In order to reproduce figure 6, it is necessary to run this script first**, creating the appropriate data structure (or adapting the script paths to the user's preference).

`05_Plots_Results.R`
Produces the plots seen in the paper. It `sources` the functions and data from the previous scripts. It also contains some specific results such as the ratio computations.

`Data_for_plots.R`
Just data to use in the plots, so that people do not have to do the whole process.

## REPRODUCIBILITY

In order to develop all the simulations (including the simulated dataset, as well as the dataset for Japan, Denmark and Iberia), it has taken us roughly a month computing in parallel with 50 cores. For a quick reproducibility check, the user must change change some lines in the code. Please note that this will not give good nor reliable results, but it can be used to make sure that the code works before launching more time consuming simulations.

In order to do this, the user must change the following, at the `03.1_ABC_Simulated.R` script:

- Line 26: Reduce the number of simulations for the object `sim`. 
- Lines 49-50: Either eliminate the loop to randomise the start or decrease strongly the number of iterations (perhaps `n_sim <- c(1:2)`?)
- Line 141: Set the desired number of cores. Personal computers will very unlikely have 50 cores.
- Line 207: Reduce the number of particles (perhaps `part_sim <- 500`?)
- Line 217: Set one single stage with while `(n_theta <= 1)`

This will only produce a ver small number of (non-usable) simulations, but this ensures that the code works and the user can then set higher values to send to a computer cluster. 

If the user wants to try two or more particles, but still run the code fast(ish), it will have to change the `line 203` and set a very high threshold. Too low unassessed thresholds might result in the stuck of the algorith, since in the parallelisation the model might never be able reach those values (hence using quantiles here).

If the user wants to do the same operation with the scripts `03.1_ABC_Observed.R` the protocol is the same, although the numbers of the lines can vary slightly. 


