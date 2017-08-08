# tropical_precip_indicators

README

This readme contains a description of the files deposited in this directory. These files were used to calculate and analyze the nine indicators introduced in the study "Quantifying the zonal-mean structure of tropical precipitation" submitted to GRL by Popp and Lutsko (2017). 

Please note that the code in this directory comes without any warranty, without even the implied warranty of merchantibility or fitness for a particular purpose. 

The directory contains four folders:
1. The folder AMIP contains the scripts used to calculate and analyze the nine indicators for 24 AMIP experiments.
2. The folder HISTORICAL contains the scripts used to calculate and analyze the nine indicators for 26 CMIP-HISTORICAL experiments. 
3. The folder AQUA contains the scripts used to calculate and analyze the nine indicators for 8 AQUA-planet experiments.
4. The folder statistical_analysis contains scripts to perform some statistical analysis of the nine indicators.

The four folders contain three types of files:
1. The folder statistical_analysis contains python scripts that perform some statistical analysis:
2. The folders AMIP, HISTORICAL and AQUA contain netcdf-files that contain the indicators calculated for the respective experiments.
3. The folders AMIP, HISTORICAL and AQUA contain ncl-scripts to perform vairous calculations related to the indicators. Note that the paths of the input files have to be adjusted for using the scripts. It is also necessary to download the simulations from the CMIP5 archives in order to perform the analysis. The scripts require temporal and zonal means of the input data.

NCL-scripts:
1. The files 'asymmetry_mlr.ncl' perform a multiple linear regression of a specified subset of the nine indicators the tropical hemispheric asymmetry index and returns the correlation between multiple linear regression and the asymmetry index. 
2. The files starting with 'cor_indicators' calculate the correlation of the indicators between different sets of experiment (for example between AMIP and CMIP-Historical).
3. The files starting with 'indicators_9' followed by a year return the correlations between the indicators for the different experiments. 
4. The files starting with 'indicators_9_mean_variance' followed by a year return the mean and the standard deviation of the nine indicators as well as the mean of the ERA-interim and GPCP data and the RMSE of the simulations compared to ERA-interim and GPCP.
5. The files starting with 'indicators_9_norm' followed by a year return the correlations between the normalized indicators for the different experiments.
6. The files starting with 'indicators_9_norm_previous' followed by a year return the correlation between the normalized indicators and previously used indicators.
7. The files starting with 'indicators_9_previous' followed by a year return the correlation between the nine indicators and the previously used indicators.
8. The files starting with 'indicators_projection_9' followed by a year or experiment types make use of equations (10) and (11) in Popp and Lutsko 2017 to approximate the zonal-mean of precipitation and creates a plot with the performance of the method.
9. The files starting with 'process_indicators' followed by an experiment type create the netcdf-files containing the indicators.

Please note that the ncl-scripts have not been cleaned up and may contain unused or redundant code. The notation of the indicators between the paper and the code here is slightly different. P_E is simply E, P_D simply D, P_A simply A, P_H simply H, P_R simply R, phi_I simply I, phi_S simply S, phi_L simply L and phi_W simply W.  
