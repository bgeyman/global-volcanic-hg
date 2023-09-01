# Global volcanic mercury (Hg) emissions based on satellite SO<sub>2</sub> measurements
## Introduction
This repository contains source code for the paper:
*Impacts of Volcanic Emissions on the Global Biogeochemical Mercury Cycle: Insights from Satellite Observations and Chemical Transport Modeling*.

Emission files used in GEOS-Chem modeling are stored separately on <a href="https://doi.org/10.7910/DVN/KHP4KK">Harvard Dataverse</a> and can be downloaded to the repository directory using *Fetch_Dataverse_Files.py*.

The script `call_all.sh` executes the python files listed below. The script `install_environment.sh` creates a conda environment called **volcano_env** that installs python and all dependencies needed to run the code.

#### Recommended citation:

> Geyman, B.M., C.P. Thackray, D.J. Jacob, and E.M. Sunderland (2023). New Satellite Data for SO<sub>2</sub> Suggest Greater Effects of Volcanic Mercury Emissions in the Northern Hemisphere. *Geophysical Research Letters*. <a href="">DOI</a>

__________

### Contents:

| Script Name | Description | Notes |
| --- | --- | --- |
| *0-Fetch_Dataverse_Files.py*  | Downloads emission files and GEOS-Chem output from Harvard Dataverse. | This is a large download (~1.6 Gb) and takes a few minutes. Files written to /Output/Publication_Data/ |
| *1-Create_Emission_Pulse_Files.py*  | Creates 100 Mg Hg emission files used in plume altitude and speciation experiment. | Currently commented out in *call_all.sh*. Uncomment to write files to /Output/Pulse_Emission_Fields/  |
| *2-Merge_SO2_Datasets.py* | Merges effusive and explosive eruption emissions from the <a href="https://disc.gsfc.nasa.gov/datasets/MSVOLSO2L4_4/summary">Multi-Satellite Volcanic Sulfur Dioxide L4 Long-Term Global Database V4</a> with passive degassing emissions from <a href="https://doi.org/10.5194/essd-15-75-2023">Fioletov et al. (2023)</a>. | Creates Output/all_SO2.csv |
| *3-Bootstrap_Hg_Emissions.py* | Generates bootstrap Hg emission estimates based on uncertainty in SO<sub>2</sub> emissions and empirical distribution function of Hg:SO<sub>2</sub> given by compiled observations. | Writes output to /Output/Data/Bootstrap_Output/ |
| *4-Estimate_Hg_Uniform.py* | Creates Hg estimates based on uniform index ratio from Bagnato et al. (2015). Reads SO<sub>2</sub> emissions from all_SO2.csv and propagates uncertainty in Hg and SO<sub>2</sub>. | |
| *5-Make_Hg_Emission_Tables.py* | Creates Table 1 in main manuscript. | Writes output to /Output/Tables/ |
| *6-Combine_RF_Hg_Emissions.py* | Adds bootstrap Hg estimates to row format inventory data | Combined output written to /Output/Data/ |
| *7-Analyze_Emission_Dist.py* | Creates figure showing meridional distribution of emissions (Fig. S1) | |
| *8-Analyze_Model_Output* | Prints budget diagnostics and creates map of Hg concentration and deposition enhancements attributable to primary volcanic emissions (Fig. 1) | Custom colormaps saved in /Data/color_lists.json |
| *9-Plot_Pulse_Experiment.py* | Plots results of plume altitude and speciation experiment (Fig. 2) | |
| *10-sevenbox_volcano_Hg.py* | Calculates steady state magnitude of atmosphere and upper ocean under volcanic Hg emissions estimated in this work. | |

__________

### Installation Instructions

### Mac
 - Install miniconda [**<a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html">Link</a>**]
 - Download and unzip `Volcanic_Hg_Cycling` repository
 - Open the terminal. At the terminal prompt, type:
   ```
   cd ~/Downloads/Volcanic_Hg_Cycling/Scripts/
   bash install_environment.sh
   conda activate volcano_env
   bash call_all.sh
   ```


