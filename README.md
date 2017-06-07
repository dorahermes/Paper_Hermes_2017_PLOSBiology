# Paper_Hermes_2017_PLOSBiology
code to download data and generate figures from Hermes, Nguyen and Winawer 2017, PLOS Biology

This repository is associated with the publication below.
Hermes, Nguyen and Winawer (2017). Neuronal synchrony and the relation between the BOLD response and the local field potential. PLOS Biology

The code downloads the data from the Open Science Framework project page (https://osf.io/v22bm/ or doi:10.17605/OSF.IO/V22BM) and then reproduces the main figures in the paper, as well as some supplementary figures.

The code runs in Matlab and was written and tested using Matlab2014b and Matlab 2017a on Apple computers with (OS 10.12.2 and 10.12.5).

Dependencies:

curve_fitting_toolbox
optimization_toolbox
signal_toolbox
statistics_toolbox
To use the code, download or clone this github repository, then navigate to the repository in the Matlab command window. Example usage below:

 % Navigate to the correct folder and add the paths
 addpath(genpath(pwd))
 
 % Download the data (this only needs to be done once)
 boldlfp_downloadData();
 
 % Make the figures
ns_script01_Fig1.m
ns_script02_Fig2.m
ns_script03_Fig3.m
ns_script05_Fig5.m

ns_script06A_RunCalibration.m 

% !!! IMPORTANT NOTES !!! 
% 1) This script takes a little while to run
% 2) This script saves ~8GB of simulated data.

ns_script06B_Fig6.m

ns_script07A_RunSimulation.m

% !!! IMPORTANT NOTES !!! 
% 1) This script takes a long time to run, and fits 8 models in 22 electrodes.
% 2) This script saves ~37GB of simulated data.

ns_script07B_Fig7AB_Fig8CD.m
ns_script07C_Fig7C_S7_S9.m
ns_script07D_Fig4.m
ns_script09A_Fig8AB_Fig9AB.m
ns_script09B_Fig9CD.m
ns_script10_Fig10.m


Code (c) Dora Hermes and Jonathan Winawer

Please direct any comments about the code via the issues page associated with this GitHub repository, or via email to dorahermes@gmail.com
