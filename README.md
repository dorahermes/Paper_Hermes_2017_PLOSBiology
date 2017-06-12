# Paper_Hermes_2017_PLOSBiology
code to download data and generate figures from Hermes, Nguyen and Winawer 2017, PLOS Biology

This repository is associated with the publication below.
Hermes, Nguyen and Winawer (2017). Neuronal synchrony and the relation between the BOLD response and the local field potential. PLOS Biology

The code downloads the data from the Open Science Framework project page (https://osf.io/v22bm/ or doi:10.17605/OSF.IO/V22BM) and then reproduces the main figures in the paper, as well as some supplementary figures.

The code runs in Matlab and was written and tested using Matlab2014b and Matlab 2017a on Apple computers with (OS 10.12.2 and 10.12.5).

Dependencies: <br/>
curve_fitting_toolbox <br/>
optimization_toolbox <br/>
signal_toolbox <br/>
statistics_toolbox <br/>

To use the code, download or clone this github repository, then navigate to the repository in the Matlab command window. Example usage below: <br/>
 % Navigate to the correct folder and add the paths <br/>
 addpath(genpath(pwd))
 
 % Download the data (this only needs to be done once) <br/>
 boldlfp_downloadData();
 
 % Make the figures <br/>
ns_script01_Fig1.m <br/>
ns_script02_Fig2.m <br/>
ns_script03_Fig3.m <br/>
ns_script05_Fig5.m <br/>

% **IMPORTANT NOTES**  <br/> 
% 1) script06A takes a little while to run (about 1 hour on a fast macbook pro)  <br/>
% 2) script06A saves ~8GB of simulated data.  <br/>

ns_script06A_RunCalibration.m  <br/>
ns_script06B_Fig6.m  <br/>

% **IMPORTANT NOTES**   <br/>
% script07A takes a long time to run (~12 hours on a fast macbook pro), and fits 8 models in 22 electrodes.  <br/>
% script07A saves ~37GB of simulated data.  <br/>

ns_script07A_RunSimulation.m  <br/>
ns_script07B_Fig7AB_Fig8CD.m  <br/>
ns_script07C_Fig7C_S7_S9.m  <br/>
ns_script07D_Fig4.m  <br/>
ns_script09A_Fig8AB_Fig9AB.m  <br/>
ns_script09B_Fig9CD.m  <br/>
ns_script10_Fig10.m  <br/>


Code (c) Dora Hermes and Jonathan Winawer

Please direct any comments about the code via the issues page associated with this GitHub repository, or via email to dorahermes@gmail.com
