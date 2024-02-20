This project models the iceflow for tidally locked exoplanets with floating ice coverage. 

The main function Floating_ice_sheets.m is modified from Tziperman et al. 2012. Additional files have supporting functions described below.

exp_00.m sets up the experiment. Additional files can be created to run multiple experiments in an array.

make_continents.m is used to create continental configurations.

output_analysis.m has some plotting functions to analyze the results output from Floating_ice_sheets.m

set_parameters.m is used to create the par structure and is called in Floating_ice_sheets.m. It shouldn't need to be edited unless new constants are added.

Parameters are saved at the end of each run from Floating_ice_sheets.m and saved under /Restart in a .mat file.
Figures generated in Floating_ice_sheets.m are saved under /Figures as PDF. The frequency of plotting at different timesteps is defined in the exp_00.m file.
