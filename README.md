# Aerosol Diagnostics Tool

**The Aerosol Diagnostics Tool produces diagnostic figures comparing simulated aerosol optical depth against remotely-sensed observational estimates. Total (AOD), absorbing (AAOD), and dust (DOD) optical depth, and angstrom exponent (AE), analyses are currently supported.** 

Three types of figures can be produced: 
  
  * **Global properties** (timeseries, seasonal cycles, zonal means, and taylor diagrams) for each variable of interest, comparing the selected models against each other and against observations. (*comparemodels_globalproperties.py*)
  * **Spatial maps** showing the distribution of each variable in a specified model, in two observational reference datasets, and the difference between these, for each season. (*singlemodel_seasonalmaps.py*)
  * **Land vs ocean properties** (timeseries and taylor diagrams) for each variable of interest, comparing the selected models against each other and against observations. (*comparemodels_landvsocn.py;* note that this feature is less polished than the others.)

The Aerosol Diagnostics Tool consists of a main coordinating script which calls separate analysis scrips, three config files, and a small number of utility scripts.  

  * The **main coordinating script** is *aerosol_diagnostics.py*. This script reads in user-specified parameters from config.yaml and outputs aerosol diagnostic figures in a series of multi-page pdfs. The wiki includes sample output pdfs for reference. It should run without specific user input, but you can edit the calls to the analysis scrips if you would like to modify what is being plotted.
  
  * The **config files** need to be filled out by the user. There are three: 
      * *config.yaml* specifies the models to be evaluated, the years + variables to be compared, types of figures to plot, and the name of the output pdf. 
      * *datapaths.yaml* specifies the paths to the different datasets.
      * *modelvars.yaml* specifies the variable names used in the simulated datasets.
  
  * The **analysis scripts** do the actual processing and plotting of the data. They are called from *aerosol_diagnostics.py* and should run without user input.
  
  * The **utility scripts** are mainly responsible for reading and processes observed and simulated optical depth fields for use in the analysis script. They should run without any user input, although *datareaders.py* does make certain assumptions about directory structure. These assumptions are outlined in the wiki. One utility script, *taylorDiagram.py*, is used for plotting. It was developed by Yannick Copin, and all rights belong to him. The original can be found at https://gist.github.com/ycopin/3342888. 
  

**High-level instructions for using the Aerosol Diagnostics Tool**
  1. download all of the observational data following the guidelines in the wiki
  2. fill out the config files
  3. execute aerosol_diagnostics.py


**The wiki provides documentation on:**
  * initial setup (required python packages, where to find the observational data, how to organize your data so the ADT can find what it needs);
  * running the ADT (filling out the config files, running the main script);
  * sample output;
  * and miscellaneous other documentation.

