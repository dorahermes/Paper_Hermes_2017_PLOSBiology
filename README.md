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
Code (c) Dora Hermes and Jonathan Winawer

Please direct any comments about the code via the issues page associated with this GitHub repository, or via email to dorahermes@gmail.com
