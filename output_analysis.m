%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis of output from Floating_ice_sheets
% and plot making
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import the variables from the output file
% this isn't for the GCM_output.nc file
FIS_output = '~/Programming/MSI-ICE-FLOW/Floating-ice-sheet-dynamics/output_history/Output-n=100000.nc';
zMOL = ncread(FIS_output,'zMol');
lon = ncread(FIS_output,'lon')';
lat = ncread(FIS_output,'lat');
time = ncread(FIS_output,'time');
%B = ncread(FIS_output,'viscosity');
%U = ncread(FIS_output,'U_veloc');
%V = ncread(FIS_output,'V_veloc');

%% Plot of maximum ice elevation
max_z = max(max(zMOL));
a = size(zMOL);
max_z = reshape(max_z, [1,a(3)]);
figure(14);clf;
% abs(max(max(z(t-1)))-max(max(z(t))))
plot(max_z(2:end),'o');
%title('','FontSize',18);
ylabel('Maximum Ice Thickness (m)','FontSize',22);
xlabel('Time (iterations)','FontSize',22);
saveas(gca, 'Convergence.png')

%% do same as above with avg Z and avg U, V in same graph
mean_Z = mean(zMOL,[1 2]);
mean_B = mean(B,[1 2]);
mean_U = mean(U,[1 2]);
mean_V = mean(V,[1 2]);
a = size(zMOL);
mean_Z = reshape(mean_Z, [1,a(3)]);
mean_B = reshape(mean_Z, [1,a(3)]);
mean_U = reshape(mean_Z, [1,a(3)]);
mean_V = reshape(mean_Z, [1,a(3)]);
figure(15);clf; % make this multiplot
plot(mean_Z(2:end),'o');

%% first derivative plot
% check difference between time steps -> convergence T-1 -T plot diff 
% then sum difference over all grid points
% save in vector and plot as a fct of time

a=size(zMOL);
der1_mat = [0 0];
for i=2:a(3)
    diff_map = zMOL(:,:,i-1) - zMOL(:,:,i);
    %diff_map = abs(diff_map);
    summed_diff = sum(diff_map,"all","omitnan");
    summed_diff = abs(summed_diff);
    der1_mat(end+1) = summed_diff;
    %der1_mat = reshape(der1_mat, [1,summed_diff]);
end

figure(16);clf;
% abs(max(max(z(t-1)))-max(max(z(t))))
plot(der1_mat(50:end),'o');
%title('','FontSize',18);
ylabel('Difference in ice elevation at each timestep (m)','FontSize',22);
xlabel('Time (iterations)','FontSize',22);
saveas(gca, 'Derivatives.png')

%% set up some parameters for the next plots
theta_north=10;
dphi=360/(176-3);
dtheta=(90-theta_north)*2/(176-3);
dphi_rad=dphi*pi/180;
dtheta_rad=dtheta*pi/180;
phi=([1:176]-2)*dphi;
theta=theta_north+([1:176]-2)*dtheta;
X = phi;
Y = 90-theta;

%% Mass balance conservation to see convergence
% verifier balance de masse (pondérée par l'aire de chaque point de grille
% create a convergence variable which calculates a criteria to assess
% convergence, calculer le volume total de glace

%% Final ice elevation
zMOL(isnan(zMOL))=0;
fig_tsurf_vectors = figure(18); clf
[~,hc]=contourf(X,Y,rot90(zMOL(:,:,end))); 
set(hc,'EdgeColor','none');
title('Ice thickness (m)','FontSize',20); 
crameri('-grayC');
xlabel('longitude','FontSize',20);
ylabel('latitude','FontSize',20);
colorbar;
saveas(gca, 'Final_ice_state_n50000.png')
colorbar;

%% make a gif of the evolution of the ice
figure(17);clf;
% 1. Create the initial image file
gifFile = 'ice_flow_50000.gif';
exportgraphics(gca, gifFile);
% 2. Within a loop, append the gif image
s = size(zMOL);
zMOL(isnan(zMOL))=0;
for i = 2:s(3)
    % Update the figure/axes
    [~,hc]=contourf(X,Y,rot90(zMOL(:,:,i)));
    set(hc,'EdgeColor','none');
    title('Ice thickness (m)','FontSize',18);
    xlabel('Longitude','FontSize',22);
    ylabel('Latitude','FontSize',22);
    crameri('-grayC');
    clim([800 3500])
    colorbar
    %annotation('rectangle',[1 1 1 1],'Color','w');
    annotation('rectangle',[0 0 1 1],'Color','w')
    exportgraphics(gca, gifFile, Append=true);
end


%% make a gif of the evolution of the U component of velocity
figure(18);clf;
% 1. Create the initial image file
gifFile = 'u_evolution.gif';
exportgraphics(gca, gifFile);
% 2. Within a loop, append the gif image
s = size(U);
zMOL(isnan(zMOL))=0;
for i = 2:s(3)
    % Update the figure/axes
    [~,hc]=contourf(X,Y,rot90(U(:,:,i)));
    set(hc,'EdgeColor','none');
    title('Ice velocity in U  (horizontal) (m/yr)','FontSize',18);
    xlabel('Longitude','FontSize',22);
    ylabel('Latitude','FontSize',22);
    crameri('vik');
    caxis([800 3500])
    colorbar
    %annotation('rectangle',[1 1 1 1],'Color','w');
    exportgraphics(gca, gifFile, Append=true);
    % change fixed axis
end

%% if needed for the GCM_output.nc file analysis
GCM_output = '~/Programming/MSI-ICE-FLOW/Floating-ice-sheet-dynamics/Input/GCM_output/surface_icedyn.001.nc';
zMOL= ncread(GCM_output,'zMOL');
tsurf = ncread(GCM_output,'tsurf');
lon = ncread(GCM_output,'longitude')';
lat = ncread(GCM_output,'latitude');

figure(1);clf;
[a,hc] = contourf(zMOL(:,:,end));
set(hc,'EdgeColor','none');
crameri('grayC')
colorbar;