%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create .dat files for continental configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grid resolution:
par.ni=176;
par.nj=176;

continent = ones(176);
continent(1:10,1:10) = 0; % set arbitrary area to 0 -> continent
save('~/Programming/MSI-ICE-FLOW/Floating-ice-sheet-dynamics/Input/masks/continent.mat','continent');
%load continent.mat   % load the file