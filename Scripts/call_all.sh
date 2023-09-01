#!/bin/bash

python 0-Fetch_Dataverse_Files.py

#python 1-Create_Emission_Pulse_Files.py

python 2-Merge_SO2_Datasets.py

python 3-Estimate_Hg_Bootstrap.py

python 4-Estimate_Hg_Uniform.py

python 5-Make_Hg_Emission_Table.py

python 6-Combine_RF_Hg_Emissions.py

python 7-Analyze_Emission_Dist.py

python 8-Analyze_Model_Output.py

python 9-Plot_Pulse_Experiment.py

python 10-sevenbox_volcano_Hg.py

# call script to remove intermediate files
python cleanup.py