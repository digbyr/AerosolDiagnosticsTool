# Aerosol Diagnostics Tool

**The Aerosol Diagnostics Tool produces diagnostic figures comparing simulated aerosol optical depth against remotely-sensed observational estimates. Total (AOD) and dust (DOD) optical depth analyses are currently supported.** 

The Aerosol Diagnostics Tool consists of a main analysis script, a data-reader utility, and a config file. 
  * The **main analysis script** is provided as both a python script (.py) and an interactive jupyter notebook (.ipynb). User input is required to specify the models, years, and variables to be plotted; the location in which these details are specified is clearly indicated in the scripts, and further explanation is provided in the wiki. Throughout the wiki, references to aerosol_diagnostics.py refer equally to the .ipynb notebook version.
  * The **data-reader utility** reads and processes observed and simulated optical depth fields for use in the analysis script. It should run without any user input, although it does make certain assumptions about directory structure. These assumptions are outlined in the wiki. All datasets, including the reference observations, need to be provided by the user.
  * The **config file** needs to be filled out by the user. It specifies the paths to the different datasets and the variable names used in the simulated datasets (see wiki for details).

**The wiki provides documentation on:**
  * The python packages required for the Aerosol Diagnostic Tool to run;
  * Instructions for filling out the config file, and expected directory structure in which data should be organized;
  * The observational datasets used (a brief overview of each instrument, and where to find the data);
  * The inputs required for each function in the data-reader and analysis scripts

