%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize all parameters for Floating_ice_sheets.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.g=10; % gravity
par.T_f=273.16; % T freeze of water
par.rho_w=1000; % water density
par.rho_I=900; % ice density
par.L_f=334e3; % latent heat of ice (J/kg),
par.mu=par.rho_I/par.rho_w; % ratio of densities
par.year=365*86400; % seconds per year
par.S0=0.015/par.year; % m/sec
par.L=2e7; % defines physical domain
par.R=6300e3; % Earth radius, assume Proxima b is earth sized, if needed

par.kappa=1.2e0;  % 1.0e0
par.T_surface_profile_type='cold';
% actual diffusivity in ice: http://www.its.caltech.edu/~atomic/snowcrystals/ice/ice.htm
% kappa=conductivity/(Cp*rho)=(2.4 [J/(sec m K)])/(1960 [J/K kg] *900 [kg/m^3])=1.36e-6 [m^2/sec] 
par.kappa_ice=1.36e-6;

par.nn=3; % exponent of Glen's law. 

% Time used to define the number of incrementations needed to span the
% value of Time
par.Time=1e5*par.year; % time span -> define nb of increments
par.dt=0.05e1*par.year; % x années en seconde
par.nt=100; % ceil(par.Time/par.dt);

par.nplot = 100; % determines frequency of plotting 
par.nsave = 10; % ceil(par.nt./par.nplot); % defines number of time steps saved as output
% par.do_use_imagesc=0; % decide between imagesc fct or contourf

par.plot.min_h=NaN; 
par.plot.max_h=NaN;
par.plot.do_h_ylabel=1;
par.nwrite_restart=min(par.nplot,floor(par.nt/10)); % freq of write restart file 

% dimensions / grid size / resolution
% options 89 / 176 -> because of Pollard forcings
par.ni=176;
par.nj=176;
par.nk=42;

% land mask initializes to 1
par.which_mask = 1; % decides between 0 land (stable) or 1 open ocean (dynamic) mask
par.land_mask=ones(par.ni,par.nj);
par.ocean_mask=ones(par.ni,par.nj);
par.no_mask=ones(par.ni,par.nj);

% following parameter chooses between a constant S taken from pollard
% (set to 1), or (if set to 0) use a constant (time-independent) S_top
% and calculate S_bot from thickness and geothermal heat flux:
par.do_constant_S=0; 

%% physical domain, including boundary points, is [2:nj-1; 2:nk-1]:
par.dzeta=1/(par.nk-3);
par.dx=par.L/(par.ni-3);
par.dy=par.L/(par.nj-3);

par.zeta=-par.mu+([1:par.nk]-2)*par.dzeta;
par.x=([1:par.ni]-2)*par.dx;
par.y=([1:par.nj]-2)*par.dy;

%% spherical coordinates:
par.theta_north=10;
par.dphi=360/(par.ni-3);
par.dtheta=(90-par.theta_north)*2/(par.nj-3);
par.dphi_rad=par.dphi*pi/180;
par.dtheta_rad=par.dtheta*pi/180;
par.phi=([1:par.ni]-2)*par.dphi;
par.theta=par.theta_north+([1:par.nj]-2)*par.dtheta;
par.s=sind(par.theta);
par.c=cosd(par.theta);

%% save in .m variable file
save('Input/parameters.mat','par');
