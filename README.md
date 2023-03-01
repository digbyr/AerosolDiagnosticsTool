# Aerosol Diagnostics Tool

**The Aerosol Diagnostics Tool produces diagnostic figures comparing simulated aerosol optical depth against remotely-sensed observational estimates. Total (AOD) and absorbing (AAOD) optical depth, and angstron exponent (AE), analyses are currently supported.** 

The Aerosol Diagnostics Tool consists of a main coordinating script which calls separate analysis scrips, three config files, and a small number of utility scripts. 

  * The **main coordinating script** is *aerosol_diagnostics.py*. This script reads in user-specified parameters from config.yaml and outputs aerosol diagnostic figures in a series of multi-page pdfs. The wiki includes sample output pdfs for reference. It should run without specific user input, but you can edit the calls to the analysis scrips if you would like to modify what is being plotted.
  
  * The **config files** need to be filled out by the user. There are three: 
      * *config.yaml* specifies the models to be evaluated, the years + variables to be compared, and the name of the output pdf. 
      * *datapaths.yaml* specifies the paths to the different datasets.
      * *modelvars.yaml* specifies the variable names used in the simulated datasets.
  
  * The **analysis scripts** do the actual processing and plotting of the data. They are called from *aerosol_diagnostics.py* and should run without user input.
  
  * The **utility scripts** are mainly responsible for reading and processes observed and simulated optical depth fields for use in the analysis script. They should run without any user input, although *datareaders.py* does make certain assumptions about directory structure. These assumptions are outlined in the wiki. One utility script, *taylorDiagram.py*, is used for plotting. It was developed by Yannick Copin, and all rights belong to him. The original can be found at https://gist.github.com/ycopin/3342888. 
  

**High-level instructions for using the Aerosol Diagnostics Tool**
  1. download all of the observational data following the guidelines in the wiki
  2. fill out the config files
  3. execute aerosol_diagnostics.py


**The wiki provides documentation on:**
  * The python packages required for the Aerosol Diagnostic Tool to run;
  * Instructions for filling out the config files, and expected directory structure in which data should be organized;
  * The observational datasets used (a brief overview of each instrument, and where to find the data);
  * The inputs used by each function in the data-reader and analysis scripts.
  * Sample output demonstrating the type of figures created by the Tool.

