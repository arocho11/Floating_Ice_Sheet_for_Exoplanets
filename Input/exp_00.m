%% Define experiment parameters
par.model_to_run='2d-sphere';
par.T_surface_profile_type='mine'; % options are mine, cold, warm
% mine is the output from icedyn (GCM), accounts for tidal locking

% import file with E-P/ melt-freeze forcing and GCM output results used as input data
% import GCM data used as input
par.S_filename='~/Programming/MSI-ICE-FLOW/Floating-ice-sheet-dynamics/Input/Pollard/pollard_forcing_interpolated_nj=176.mat';
par.GCM_input = '~/Programming/MSI-ICE-FLOW/Floating-ice-sheet-dynamics/Input/GCM_output/surface_icedyn.001.nc';

% (to be edited) option to name the output file 
% par.GCM_filename = '';

par.dt=0.05e1*par.year; % x years to seconds
par.nt=10; 

par.nplot = 75; % determines frequency of plotting (plots every nplot)
par.nsave = 2; % determines frequency of saving in zMol output mat

%% mask:

% determine which type of mask to use
% options are 0 -> static continental land_mask; 1 -> dynamic ocean_mask; 2 -> both
par.which_mask = 2;

% grid resolution:
par.ni=176;
par.nj=176;

% if land_mask is needed, use following code
if par.which_mask == 0 || par.which_mask == 2
    % either use the continental maps given by pollard or .mat file
    fid=fopen('~/Programming/MSI-ICE-FLOW/Floating-ice-sheet-dynamics/Input/masks/630Ma_mask_x2.dat','r');
    mask=fscanf(fid,'%1d',[par.ni,par.nj]);
    par.land_mask=mask(:,:);
    fclose('all');

    % or use a personalized continental mask created using make_continents.m
    % uncomment following lines to do so:
    % load Input/masks/continent.mat
    % par.land_mask = continent;

    % it can also be set to cover open ocean from GCM input
    % just comment this whole section
end