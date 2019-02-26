options(error = utils::dump.frames)

set.seed(0)

#################################################################################
#### packages used
#################################################################################
library(parallel)
library(BASS)
library(mnormt)
library(ncdf4)
library(ggmap)
library(ggplot2)
library(RgoogleMaps)
library(fields)
library(wesanderson) # available using devtools::install_github("karthik/wesanderson")
library(MASS)
library(dismo)
library(fdakma)
library(tgp)
library(BART)
library(randomForest)
library(truncnorm)

#################################################################################
#### number of cores (or threads)
# Note: more cores provides substantial speedup in emulation, sensitivity analysis, 
# and emulator comparison, where many tasks are embarassingly parallel.  The 
# calibration portions also benefit from multiple cores, though the speedup 
# diminishes after about 4 cores.
#################################################################################
ncores<-detectCores() # number of cores


#################################################################################
#### Read in Data
#################################################################################
t1<-proc.time()
source('code/makeData_sim.R') # get simulation data, including EOFs
source('code/makeData_obs.R') # get observation data
(data.time<-proc.time()-t1)

# note: for a smaller memory footprint, we delete after each step (important quantities are saved to file)
rm(list=ls()); gc(); ncores<-detectCores()

#################################################################################
#### Build emulator
#################################################################################
t1<-proc.time()
source('code/emulator/eofBASS.R') # fit BMARS emulator in EOF space (1 BMARS model for each EOF)
(emu.time<-proc.time()-t1)
source('code/emulator/emulatorPlots.R') # plot some results

rm(list=ls()); gc(); ncores<-detectCores()

#################################################################################
#### Emulator sensitivity
#################################################################################
t1<-proc.time()
source('code/emulator/sensitivity.R') # Sobol sensitivity analysis of emulator (functional)
(sens.time<-proc.time()-t1)

rm(list=ls()); gc(); ncores<-detectCores()

#################################################################################
#### Calibration to simulated data
#################################################################################
t1<-proc.time()
source('code/calibration/calibSim.R') # calibrate to one of the holdout model runs
(calibSim.time<-proc.time()-t1)
source('code/calibration/calibPlots.R') # plot results

#################################################################################
#### Calibration to actual data
#################################################################################
t1<-proc.time()
source('code/calibration/calibObs.R') # calibrate using observations
(calib.time<-proc.time()-t1)
source('code/calibration/calibPlots.R') # plot results
source('code/calibration/discrepPlots.R') # plot discrepancy

rm(list=ls()); gc(); ncores<-detectCores()

#################################################################################
#### Comparing emulators
#################################################################################
t1<-proc.time()
source('code/emulator/compare/tgp_bart_bass.R') # friedman boolean simulation
source('code/emulator/compare/tgp_bart_bass2.R') # subsampling real dataset simulation
(sim.time<-proc.time()-t1)

