# Inferring-Atmospheric-Release-Characteristics-in-a-Large-Computer-Experiment-using-Bayesian-Adaptive

# Author Contributions Checklist Form

## Data

### Abstract 

An ensemble of FLEXPART-WRF particle dispersion model runs with input settings corresponding to different hypothetical releases near the Diablo Canyon Nuclear Power Plant. Measurements of tracer concentrations of an experimental release.

### Availability 

Unrestricted.

### Description

The measurement data are a product of the Pacific Gas and Electric Company, and are released with permission.

The simulation data were created under work funded by Laboratory Directed Research and Development projects at the Lawrence Livermore National Laboratory. The work was performed under the auspices of the US Department of Energy by Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344, and is released under UCRL number LLNL-MI-733093.

The data will be provided via dataverse.org.

The simulation data set contains an ensemble of 18000 runs of FLEXPART-WRF, a particle dispersion model. The model output are time series of sulfur hexafluoride concentrations at 137 locations near the Diablo Canyon Nuclear Power Plant. The model takes 11 inputs, six corresponding to FLEXPART inputs and five corresponding to WRF inputs.

The measurement dataset contains time series of sulfur hexafluoride concentrations at 150 locations after an experimental release from Diablo Canyon Nuclear Power Plant in 1986.

Using these data, the goals are to (1) use the simulations to construct a statistical emulator of the particle dispersion model and (2) use the measurements and the emulator to identify the most likely inputs to the simulator.

The simulated data are contained in the file flxout-dopptex-timeseries-lhs04.nc, with the inputs to the simulations given in flexuq_lhs04-wrfdims.tab. The measurement data are contained in the DOPPTEX folder, with descriptions of the data therein.

## Code

### Abstract

The code reads the simulation data, builds an emulator, performs a sensitivity analysis, reads the calibration data, and performs the calibration.

### Description 

The code will be delivered via a github repository of R code, and is is protected under the GNU General Public License, version 3 (GPL-3). The code relies heavily on the R package BASS, version 0.2.2, and less heavily on a number of other packages. The code:
* reads the simulation and observation data in their original formats;
* gets the EOF decomposition of the simulation data (100 EOFs used);
* builds 100 BMARS models for the EOF weights;
* performs a functional sensitivity analysis of the emulator;
* performs a synthetic calibration to one of the model runs;
* performs calibration to field measurements with discrepancy; 
* performs a functional clustering analysis of the discrepancy;
* does an emulator comparison simulation.

## Instructions for Use

### Reproducibility 

To reproduce figures 3-7 and 9-11 in the manuscript, as well as the three figures in the supplement, use the script wrapper.R. While the code to reproduce figures 1, 2, and 8 is also included, it requires Google Maps API access. The exact results are random seed dependent, but without a matching seed the same general results are obtained.

Some of the pieces of the code used for this analysis make extensive use of parallelism and have moderately large memory footprint (especially when using parallelism). To achieve some sense of modularity, the different tasks of the code (such as reading in the data, building the emulator, performing calibration, etc.) are in separate files that use the output saved from any relevant upstream processes (for instance, the calibration portion uses the emulation output, but not the sensitivity analysis). Following is a list of the parts of the code and the time they take with the preferred computational resources.

* **Reading/processing data:** 25 minutes (one core, but faster with multithreaded BLAS).
* **Building emulator:** the longest running models take about two hours, and the 100 models can be fit in parallel, so there is opportunity for this to take only as long as it takes to fit the longest models. We use 64 cores, and that is enough so that all models are fit by the time the longest one finishes (i.e., jobs are not prescheduled).
* **Emulator comparison:** Using the same number of cores as simulations, this takes about 20 minutes.

### Additional Notes

Since the code and data are not yet publicly available, use the zipped file included in the submission.
