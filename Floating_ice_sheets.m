function Floating_ice_sheets(expnum)
% Solve the model for Marine ice sheets adapted from Tziperman et al. 2012, Eli, June 2011
% for exoplanets: tidally-locked floating ice covered water world, expected to have 
% an open ocean at the substellar point
% Alexandra Rochon, undergraduate researcher 2023

% experiment parameters are defined in a separate exp_00.m file conventions:
% i,j is longitude, latitude
% u pos towards the east, (west-> east)
% v pos towards equator (from north goes down / from south goes up)
 
par=set_parameters(expnum);
integrate_h_2d_sphere(par);

fprintf(1,'done.\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par=set_parameters(expnum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define base parameters, overwrite with params from exp_XX.m file where
% needed, some standard parameters are defined in set_parameters.m and
% saved in paramaters.mat variable file

addpath Input;  % Could also put this line in my ~/matlab/startup.m
addpath Restart;
load('Input/parameters.mat', 'par')

par.expnum=expnum;

% read scenario parameters:
% takes in parameters from exp_00.m file
eval(sprintf('exp_%2.2d',par.expnum));

%% Import values zMol and Surface Temp from GCM results
% this could be commented out or taken into exp_03 if the data files cant
% be found, this is used to couple with GCM
zMOL = ncread(par.GCM_input,'zMOL');
tsurf = ncread(par.GCM_input,'tsurf');
lon = ncread(par.GCM_input,'longitude')';
lat = ncread(par.GCM_input,'latitude');

X = linspace(lon(1),lon(end),176);
Y = linspace(lat(end),lat(1),176)';

% interpolate GCM data to [176, 176] resolution and get all values w.r.t.
% 0 instead of relative elevation
% axes are reversed to follow the (lon,lat) convention established here
newZ = interp2(lat,lon,zMOL,Y,X);
par.minZ = min(min(newZ));
newZ = (newZ.*1000)+abs(par.minZ.*1000);
% create mask covering any initial ice thickness below 10 m
% par.land_mask(newZ<10) = 0; % to use static land_mask covering open ocean 
par.ocean_mask(newZ<10) = 0;
par.zMol = newZ; 
% interpolate surface temp
newTS = interp2(lat,lon,tsurf,Y,X);
par.tSurf = newTS;
% average surface temp in both dir
par.tLat = mean(newTS, 1)'; % used in plotting
par.tLon = mean(newTS, 2)';
par.B_ocean = 1e16.*mean(par.zMol,'all')./par.R; % average viscosity value


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,u_n,v_n,B]=read_restart(par,which_restart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restart_filename=sprintf('Restart/restart-exp-%2.2d-%s.mat',par.expnum,which_restart);
if exist(restart_filename,'file')
  load(restart_filename,'h','u_n','v_n','mask','B');
  par.ocean_mask = mask; % adjust mask here ?

else  
  fprintf(1,'no restart file %s, initializing with defaults.\n',restart_filename);
  u_n=ones(par.ni,par.nj)/par.year; % initializes velocities to 
  v_n=ones(par.ni,par.nj)/par.year; % 1 m / seconds per year
  h = zeros(par.ni,par.nj,3)+par.zMol; % initial zMol value
  B=ones(par.ni,par.nj).*1e16.*par.zMol./par.R; % newton * s / m^2
  
  B(par.ocean_mask==0) = par.B_ocean; % make sure no 0 elements in B
  B(B==0) = par.B_ocean; % set to an average value, could be a min value instead
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function integrate_h_2d_sphere(par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ice flow in two horizontal dimensions using spherical coordinates.

fprintf(1,'running nonlinear 2D-sphere model, exp num %2.2d\n',par.expnum);
fprintf(1,'Number of timesteps: %2.2d\n',par.nt);
fprintf(1,'Time incremment between each steps: %4.2f years \n',par.dt/par.year);

if 0
  % print diagnostics regarding the time step size:
  fprintf(1,'(par.R*par.s(par.nj-1)*(par.dphi_rad))^2/kappa=%g, dx/(U=50m/yr)=%g, dt=%g (yrs)\n' ...
          ,((par.R*par.s(par.nj-1)*(par.dphi_rad))^2/par.kappa)/par.year ...
          ,((par.R*par.s(par.nj-1)*(par.dphi_rad))/(50/par.year))/par.year ...
          ,par.dt/par.year);
end

[h,u,v,B]=read_restart(par,'2d-sphere-nonlinear');

%% melting/ freezing forcing:
% S=S_total(par); % defined by Pollard calcs. set to zero for now
% longitunal avg so not relevant for this case
% make sure that zonally integrated forcing over the domain is equal to S even when there are land masses:
S=zeros(par.ni,par.nj);

% uncomment this section if S value is no longer 0
% for j=1:par.nj
%   S(:,j)=ones(par.ni,1)*S(j);
% end
% mean_S=0;
% weight=0;
% for j=2:par.nj-1
%   for i=2:par.ni-1
%   mean_S=mean_S+S(i,j)*par.s(j)*par.land_mask(i,j);
%   weight=weight+par.s(j)*par.land_mask(i,j);
%   end
% end
% S(1,:)=S(par.ni-1,:);
% S(par.ni,:)=S(2,:);
% mean_S=mean_S/weight;
% S=(S-mean_S).*par.land_mask;

%% initialize variables:
div_hv_n=zeros(par.ni,par.nj);
div_hv_np1=zeros(par.ni,par.nj);
kappa_del2_h_n=zeros(par.ni,par.nj);
kappa_del2_h_np1=zeros(par.ni,par.nj);
u_np1=u;
v_np1=v;

% save every nsave iteration all values, last layer could be
% empty if neg h is found
mat_Z = zeros(par.ni,par.nj); zeros(par.ni,par.nj); 
mat_U = zeros(par.ni,par.nj);
mat_V = zeros(par.ni,par.nj);
mat_B = zeros(par.ni,par.nj);

% par.land_mask = continental configuration (static)
% par.ocean_mask = open ocean (dynamic)
if par.which_mask == 0
    mask=par.land_mask;
elseif par.which_mask == 1
    par = adjust_mask(par, h(:,:,3));
    mask=par.ocean_mask;
elseif par.which_mask == 2  
    par = adjust_mask(par, h(:,:,3));
    mask = par.land_mask.*par.ocean_mask;
else
    par.no_mask = ones(par.ni,par.nj); 
    mask=par.no_mask;
end

%% Time step h-equation:
% number of iterations to change to acceleration
for n=1:par.nt
  par.n=n; % store incrementation variable in par struct to use in fct calls
  time_kyr=n*par.dt/par.year/1000; 
  % every incrementation increases by dt, n is nb of incrementations
  v_n=v_np1;
  u_n=u_np1;
  T_surface=T_surface_function(par);
  
  %% calculate velocity:
  max_iter_eff_viscosity=20;
  par.iter_eff_viscosity=0;
  diff_eff_viscosity=1;
  count = 1;
  % viscosity calculation and convergence loop
  while par.iter_eff_viscosity<=max_iter_eff_viscosity && diff_eff_viscosity>1.e-4
    % Average viscosity over relevant temperature range:
    B_old=B;
    B=calc_eff_viscosity(par,T_surface,h(:,:,2),u_np1,v_np1);
    [u_np1,v_np1]=solve_uv_2d_matrix_form_sphere(par,B,squeeze(h(:,:,2)));
    
    %%% CHECK HERE WHAT THE IMPACT OF WHICH MASK HAS WITH MINIMAL AVG B VALUE
    
    if par.which_mask == 0
        max_B=max(max(B(find(par.land_mask~=0))));
    elseif par.which_mask == 1
        max_B=max(max(B(find(par.ocean_mask~=0))));
    elseif par.which_mask == 2
        total_mask = par.ocean_mask.*par.land_mask;
        max_B=max(max(B(find(total_mask~=0))));
    else
        max_B=max(max(B));
    end

    diff_eff_viscosity=max(max(abs(B-B_old)))./max_B;
    par.iter_eff_viscosity=par.iter_eff_viscosity+1;
    count = count+1;
    if par.iter_eff_viscosity==max_iter_eff_viscosity
      fprintf(1,'*** at n=%d, reached max_iter_eff_viscosity=%d, with diff_eff_viscosity=%g\n',n,max_iter_eff_viscosity,diff_eff_viscosity);
    end
  end

  %% set v=0 at the southern & northern boundaries (pole):   
  u_np1(:,2)=0; v_np1(:,2)=0;
  u_np1(:,par.nj-1)=0; v_np1(:,par.nj-1)=0;

  if n==1 
    u_n=u_np1; v_n=v_np1; 
  end

  %% calculate S_bot for this time step:
  if par.do_constant_S==1
    S_bot_n=0*h(:,:,1);
    S_bot_np1=0*h(:,:,2);
  else % if pollard S isn't used, modify here to account for 
      % geothermal heat flux, for now set to 0
    S_bot_n=S_bot(par,h(:,:,1));
    S_bot_np1=S_bot(par,h(:,:,2));
  end
  
  %% advance h in time using second order Adams Bashforth:
  for j=2:(par.nj-1)
    s_jmhalf=0.5*(par.s(j)+par.s(j-1));
    s_jphalf=0.5*(par.s(j)+par.s(j+1));
    for i=2:(par.ni-1)
      kappa_del2_h_n(i,j)=par.kappa*( ...
        (1/(par.R^2*par.s(j)))* ...
        ((h(i+1,j,1)-h(i,j,1))*mask(i,j)*mask(i+1,j) ...
         -(h(i,j,1)-h(i-1,j,1))*mask(i,j)*mask(i-1,j)  ...
         )/par.dphi_rad^2+...
        (1/(par.R^2*par.s(j)))* ...
        (s_jphalf*(h(i,j+1,1)-h(i,j,1))*mask(i,j)*mask(i,j+1) ...
         -s_jmhalf*(h(i,j,1)-h(i,j-1,1))*mask(i,j)*mask(i,j-1) ...
         )/par.dtheta_rad^2);
      
      div_hv_n(i,j)= (1/(par.R*par.s(j)))*(...
          (h(i+1,j,1)*u_n(i+1,j)-h(i-1,j,1)*u_n(i-1,j))/(2*par.dphi_rad)+...
          (par.s(j+1)*h(i,j+1,1)*v_n(i,j+1)-par.s(j-1)*h(i,j-1,1)*v_n(i,j-1))/(2*par.dtheta_rad));
      RHS_n=-div_hv_n(i,j)+kappa_del2_h_n(i,j)+S(i,j)+S_bot_n(i,j);

      kappa_del2_h_np1(i,j)=par.kappa*( ...
        (1/(par.R^2*par.s(j)))* ...
        ((h(i+1,j,2)-h(i,j,2))*mask(i,j)*mask(i+1,j) ...
         -(h(i,j,2)-h(i-1,j,2))*mask(i,j)*mask(i-1,j)  ...
         )/par.dphi_rad^2+...
        (1/(par.R^2*par.s(j)))* ...
        (s_jphalf*(h(i,j+1,2)-h(i,j,2))*mask(i,j)*mask(i,j+1) ...
         -s_jmhalf*(h(i,j,2)-h(i,j-1,2))*mask(i,j)*mask(i,j-1) ...
         )/par.dtheta_rad^2);

      div_hv_np1(i,j)= (1/(par.R*par.s(j)))*(...
          (h(i+1,j,2)*u_np1(i+1,j)-h(i-1,j,2)*u_np1(i-1,j))/(2*par.dphi_rad)+...
          (par.s(j+1)*h(i,j+1,2)*v_np1(i,j+1)-par.s(j-1)*h(i,j-1,2)*v_np1(i,j-1))/(2*par.dtheta_rad) ...
        );
      RHS_np1=-div_hv_np1(i,j)+kappa_del2_h_np1(i,j)+S(i,j)+S_bot_np1(i,j);

      %%%%%% check this
      h(i,j,3)=h(i,j,2)+par.dt*(1.5*RHS_np1-0.5*RHS_n); 
      
    end
  end

  %% periodic boundary conditions in x:
  h(1,:,3)=h(par.ni-1,:,3);
  h(par.ni,:,3)=h(2,:,3);
  
  %% set north and south boundary conditions of h_y=0:
  h(:,1,3)=h(:,2,3);
  h(:,par.nj,3)=h(:,par.nj-1,3);

  %% prepare for next time step:
  h(:,:,1)=h(:,:,2);
  h(:,:,2)=h(:,:,3);

  %% check for non-physical values:
  % if h goes below 0, there is a problem that the code cant handle usually
  % has something to do with velocities being so fast or timestep too large
  % that the ice flows too much, check physics to make sure the code is
  % handling everything correctly
  
  found_negative_h=0;
  %min(h(:,:,3),[],"all")
  if min(min(h(:,:,3)))<0
    fprintf(1,'*** at n=%d, min(h)=%g; max(h)=%g; \n',n,min(min(h(:,:,3))),max(max(h(:,:,3))));
    disp('*** stopping due to h<0 ***');
    found_negative_h=1;
    %min(h,[],"all")
    disp('problem!!');
  end

  % if there is any h neg it is replaced under mask
  %% should this be h(:,:,3) or (:,:,2)
  h(:,:,3)=h(:,:,3).*mask;

  %% update mask here  
  par.ocean_mask = mask;

  if par.which_mask == 0
      mask=par.land_mask;
  elseif par.which_mask == 1
      par = adjust_mask(par,h(:,:,3));
      mask=par.ocean_mask;
  elseif par.which_mask == 2
      par = adjust_mask(par,h(:,:,3));
      mask=par.ocean_mask.*par.land_mask;
  else
      mask=par.no_mask;
  end
  
%% save values in matrix
  ib=2; ie=par.ni-1; jb=2; je=par.nj-1;
  is=floor(par.ni/40); js=floor(par.nj/40);
  hplot=h(:,:,3); 
  hplot(mask==0)=NaN;
  X = par.phi(ib:ie);
  Y = 90-par.theta(je:-1:jb);
  Z = hplot(ib:ie,jb:je);
  
  u_arr = u_np1*par.year;
  v_arr = v_np1*par.year;

  if n==1 
      mat_Z = hplot;
      mat_U = u_arr; 
      mat_V = v_arr;
      mat_B = B;
  elseif floor(n/par.nsave)*par.nsave==n || n==par.nsave
      mat_Z = cat(3, mat_Z, hplot);
      mat_U = cat(3, mat_U, u_arr );
      mat_V = cat(3, mat_V, v_arr );
      mat_B = cat(3, mat_B, B );
  end 
  % progression count
  fprintf(1,'*** reached n=%d \n',n);


  if n==1 || floor(n/par.nplot)*par.nplot==n|| n==par.nt || found_negative_h==1
  %% plot all recurring figures 
  % this defines the recurrence of plotting, and plots first / last / error
  % as well
  
    if par.which_mask == 0
        mask=par.land_mask;
    elseif par.which_mask == 1
        par = adjust_mask(par,hplot);
        mask=par.ocean_mask;
    elseif par.which_mask == 2
        par = adjust_mask(par,hplot);
        mask=par.ocean_mask.*par.land_mask;
    else
        mask=par.no_mask;
    end

    %% plot (u,v) arrows and contours of thickness h:
    fig_h_vectors=figure(1); clf
    
    nan_mask=mask; 
    nan_mask(mask==0)=NaN;
    
    % c, hc returns the contour matrix and the contour object hc. Use hc to set properties
    % contourf(X,Y,Z) specifies the x and y coordinates for the values in Z.
    % X = par.phi(ib:ie)
    % Y = 90-par.theta(je:-1:jb) where par.theta is 1 x 176 double
    % Z =  hplot(ib:ie,je:-1:jb)' where hplot=h(:,:,3)
    % ib = 2, ie=par.ni-1, is=floor(par.ni/40), (j values defines same)

    [~,hc]=contourf(X,Y,Z); 
    set(hc,'EdgeColor','none');
    % crameri('-vik') % if crameri is downloaded use this colormap
    colorbar

    % to set a limit on the colorplot elevation
    if ~isnan(par.plot.min_h)
      set(gca,'clim',[par.plot.min_h, par.plot.max_h]); 
    end

    hold on
    % add velocity vector arrows with quiver plot
    quiver_scale=0.7;
    my_scale=0.005;
    sc = 2;
    
    [phi,theta]=meshgrid(par.phi(ib:is:ie),90-par.theta(je:-js:jb));
    phi=phi.*nan_mask(ib:is:ie,jb:is:je);
    x_quiver = phi(1:sc:end,1:sc:end);
    y_quiver = theta(1:sc:end,1:sc:end);
    u_quiver = my_scale*u_np1(ib:is*sc:ie,jb:is*sc:je)*par.year;
    v_quiver = my_scale*v_np1(ib:is*sc:ie,jb:is*sc:je)*par.year;
    quiver(x_quiver, y_quiver, v_quiver, u_quiver, quiver_scale,'color','k');
        
    uvmax=par.year*max(max(max(sqrt(u_np1.^2+v_np1.^2))));
    h3=title(sprintf('h;max velocity=%3.0fm/yr',uvmax));
    h1=xlabel('longitude');
    if par.plot.do_h_ylabel==1; h2=ylabel('latitude'); end

    % make sure plot size is exactly the same with and without
    % latitude label:
    set(gca,'PlotBoxAspectRatioMode','manual');
    set(gca,'Position',[0.13, 0.11, 0.718, 0.815]);

    xlim([0,360]);
    ylim([-80,80]);
    if par.plot.do_h_ylabel==1; set(h2,'fontsize',10); end
    set([gca,h1,h3],'fontsize',10);
    
    %% plot terms in h equation:
    fig_terms=figure(2); clf
    PLOT=zeros(par.ni,par.nj);
    for j=jb:je; PLOT(ib:ie,j)=mask(ib:ie,j).*S(ib:ie,j).*par.s(j); end
    hl1=plot(90-par.theta(je:-1:jb),mean(PLOT(ib:ie,je:-1:jb),1)*par.year,'r');
    hold on
    for j=jb:je; PLOT(ib:ie,j)=mask(ib:ie,j).*div_hv_n(ib:ie,j).*par.s(j); end
    hl2=plot(90-par.theta(je:-1:jb),mean(PLOT(ib:ie,je:-1:jb),1)*par.year,'--g');
    for j=jb:je; PLOT(ib:ie,j)=mask(ib:ie,j).*kappa_del2_h_n(ib:ie,j).*par.s(j); end
    hl3=plot(90-par.theta(je:-1:jb),mean(PLOT(ib:ie,je:-1:jb),1)*par.year,'--c');
    rhs=div_hv_n-kappa_del2_h_n;
    for j=jb:je; PLOT(ib:ie,j)=mask(ib:ie,j).*rhs(ib:ie,j).*par.s(j); end
    hl4=plot(90-par.theta(je:-1:jb),mean(PLOT(ib:ie,je:-1:jb),1)*par.year,'--b');
    xlim([90-par.theta(je),90-par.theta(jb)]);
    h1=xlabel('latitude');
    h2=ylabel('m/yr');
    h3=title(sprintf('terms in h eqn, n=%d, t=%gkyr',n,time_kyr));
    h4=legend('S','\nabla vh','k\nabla^2h','rhs','location','southeast');
    set([hl1,hl2,hl3,hl4],'linewidth',2)
    set([gca,h1,h2,h3,h4],'fontsize',10);

    
    %% plot effective viscosity:
    fig_B=figure(3); clf; 
    Bplot=B*par.R./h(:,:,3); 
    Bplot(mask==0)=NaN; 
    [~,hc]=contourf(par.phi(ib:ie),90-par.theta(jb:je),rot90(log10(Bplot(ib:ie,jb:je)))); set(hc,'EdgeColor','none');
    set(gca, 'YDir','reverse'); % create a mirror image of the plot to get right dir
    colorbar; 
    h1=xlabel('longitude');
    if par.plot.do_h_ylabel==1; h2=ylabel('latitude'); end
    h3=title('log_{10}(B)');
    set([gca,h1,h3],'fontsize',10);
    if par.plot.do_h_ylabel==1; set(h2,'fontsize',10); end

    
    %% plot Peclet number:
    Peclet=NaN(par.ni,par.nj);
    TS=T_surface_function(par);
    
    for j=2:par.nj-1
      for i=2:par.ni-1
            Peclet(i,j)=(v_np1(i,j)*(par.R*par.s(j))^(-1)*(par.s(j+1)*TS(i,j+1)-par.s(j-1)*TS(i,j-1))/(2*par.dtheta_rad))/(par.kappa_ice*(TS(i,j)-par.T_f)/h(i,j,3)^2);
      end
    end

    fig_Peclet=figure(4); clf; 
    Peclet(mask==0)=NaN; 
    [~,hc]=contourf(par.phi(ib:ie),90-par.theta(je:-1:jb),rot90(abs(Peclet(ib:ie,jb:je))));
    set(hc,'EdgeColor','none');
    colorbar; 
    h1=xlabel('longitude');
    h2=ylabel('latitude');
    h3=title('Peclet number');
    set([gca,h1,h2,h3],'fontsize',10);

    %% ice thickness without vectors
    fig_h_ice=figure(5); clf
    new_z = Z;
    new_z(isnan(new_z))=0;
    [~,hc]=contourf(X,Y,rot90(new_z)); set(hc,'EdgeColor','none');
    title('ice thickness (m)');
    xlabel('longitude','FontSize',20);
    ylabel('latitude','FontSize',20);
    % crameri('-grayC');
    colorbar;
    
    %% h evolution at specific location
    % doesnt produce an output at n==1 since this shows evolution
    if n~=1
        fig_h_evolution=figure(6); clf
        sz = size(mat_Z,3);
        avg_arr = zeros(1,sz);
        % take avg for all latitude at a given longitude, at each time step
        for i=1:sz
            avg_arr(1,i) = mean2(mat_Z(:,20,i));
        end
        time_arr = linspace(0, par.dt*n/par.year, sz);
        plot(time_arr, avg_arr, '--b')
        title(sprintf('h evolution, n=%d, t=%gkyr',n,time_kyr));
        xlabel('Time (years)');
        ylabel('Height of Ice (m)');

        savefig_pdf(par, fig_h_evolution, 'h-evolution');
    end

    %% u and v velocity fields
    fig_u_velocity = figure(7); clf
    [~,hc]=contourf(X,Y,rot90(u_arr(ib:ie,jb:je))); set(hc,'EdgeColor','none');
    title('horizontal velocity (u) m/yr'); 
    xlabel('longitude');
    ylabel('latitude');
    colorbar;
    
    fig_v_velocity = figure(8); clf
    [~,hc]=contourf(X,Y,rot90(v_arr(ib:ie,jb:je))); set(hc,'EdgeColor','none');
    title('vertical velocity (v) m/yr');
    xlabel('longitude');
    ylabel('latitude');
    colorbar;
        
    %% save figures as pdf using savefig_pdf function
    savefig_pdf(par, fig_h_vectors, sprintf('h-vectors_n=%2.2d',n))
    savefig_pdf(par, fig_terms, sprintf('terms_n=%2.2d',n));
    savefig_pdf(par, fig_B, sprintf('B_n=%2.2d',n));
    savefig_pdf(par, fig_Peclet, sprintf('Peclet_n=%2.2d',n));
    savefig_pdf(par, fig_h_ice, sprintf('h-ice_n=%2.2d',n));
    savefig_pdf(par, fig_u_velocity, sprintf('u-velocity_n=%2.2d',n));
    savefig_pdf(par, fig_v_velocity, sprintf('v-velocity_n=%2.2d',n));
    
    %pause
    if found_negative_h==1
      return
    end
   end % plotting

   %% Write restart
   if n==par.nt || floor(n/par.nwrite_restart)*par.nwrite_restart==n
       save(sprintf('Restart/restart-exp-%2.2d-2d-sphere-nonlinear.mat',par.expnum),'h','u_n','v_n','mask','B');
   end

  end 
    
  %% Surface temperature map with ice velocity vectors 
  % only plot once
  fig_tsurf_vectors = figure(9); clf
  [~,hc]=contourf(X,Y,rot90(TS(ib:ie,jb:je))); set(hc,'EdgeColor','none');
  title('Surface temperature (K)'); 
  % hold on
  % hq=quiver(x_quiver, y_quiver, v_quiver, u_quiver, quiver_scale,'color','k');
  % title('Surface temperature and velocity vectors (m/yr)'); 
  colormap('jet');
  xlabel('longitude','FontSize',20);
  ylabel('latitude','FontSize',20);
  colorbar;
  savefig_pdf(par, fig_tsurf_vectors, sprintf('tsurf-vectors_n=%2.2d',n));

  %% Model output and CGM input file
  % GCM input file must be final state (time step) only, but export all saved 
  % timesteps included in separate file, variables vary with time only in
  % Output.nc file, done after final iterations of the timestep loop
  
  % create all variable fields
  s = size(mat_Z);
  nccreate('Output.nc','lon','Dimensions',{'lon' 176});
  nccreate('Output.nc','lat','Dimensions',{'lat' 176});
  nccreate('Output.nc','time','Dimensions',{'time' s(3)});
  nccreate('Output.nc','zMol','Dimensions',{'lon' 176 'lat' 176 'time' s(3)});
  nccreate('Output.nc','U_veloc','Dimensions',{'lon' 176 'lat' 176 'time' s(3)});
  nccreate('Output.nc','V_veloc','Dimensions',{'lon' 176 'lat' 176 'time' s(3)});
  nccreate('Output.nc','viscosity','Dimensions',{'lon' 176 'lat' 176 'time' s(3)});

  ncwrite('Output.nc','zMol',mat_Z);
  ncwrite('Output.nc','U_veloc',mat_U);
  ncwrite('Output.nc','V_veloc',mat_V);
  ncwrite('Output.nc','viscosity',mat_B);

  % create all variable fields
  % doesnt include albedo and thermal...
  nccreate('GCM_input.nc','lon','Dimensions',{'lon' 360});
  nccreate('GCM_input.nc','lat','Dimensions',{'lat' 180});
  nccreate('GCM_input.nc','tsurf','Dimensions',{'lon' 360 'lat' 180}); % doesnt vary with time in model run
  nccreate('GCM_input.nc','zMOL','Dimensions',{'lon' 360 'lat' 180});
  
  % interpolated to GCM dim
  Xq = linspace(1,360, 360);
  Yq = linspace(1,180, 180)';
  Y_174 = linspace(1,180, 174)'; % redefine bc the Y used to plot goes from (-80,80)
  Y_176 = linspace(1,180, 176)';
  zMol = (Z/1000) - par.minZ; % might need a way to handle NaN here
  zMol = interp2(Y_174, X, zMol,Yq,Xq);
  tsurf = interp2(Y_176, par.phi, TS,Yq,Xq);

  ncwrite('GCM_input.nc','tsurf',tsurf);
  ncwrite('GCM_input.nc','zMOL',zMol);
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,v]=solve_uv_2d_matrix_form_sphere(par,B,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve diagnostic momentum equations for 2d linearized model in
% spherical coordinates, by writing them in matrix form.
% the persistent variable remain in memory even as the program stops
% runnning, s.t. every time Floating_ice_sheets(expnum) runs, it has the
% latest value of the persistent vars already stored, unless ''clear all'
% is used, in which case, all variables are initialized to '[]'
% corresponds to eqn 4, 5 in the paper
persistent im jm iuv A RHS neqns nzmax 

if par.which_mask == 0
    mask=par.land_mask;
elseif par.which_mask == 1
    par = adjust_mask(par,h);
    mask=par.ocean_mask;
elseif par.which_mask == 2
    par = adjust_mask(par,h);
    mask=par.ocean_mask.*par.land_mask;
else
    mask=par.no_mask;
end

if par.n==1 
  % define grid for first iteration
  % order variables in a single vector, including all ocean grid points
  % for both (u,v):
  neqns_max=par.ni*par.nj*2; % 1x1 double, total area
  im=zeros(neqns_max,1); jm=im; iuv=im; % initialize 3 vectors to area length to 0
end

N=zeros(par.ni,par.nj,2); % full grid, two layers
n=0;
for i=2:par.ni-1
  for j=3:par.nj-2
      if mask(i,j)==1 % skip over mask  
        % indices for u-velocities:
        n=n+1; 
        im(n)=i; jm(n)=j; iuv(n)=1;
        N(i,j,1)=n; % stores indices
        % indices for v-velocities:
        n=n+1;
        im(n)=i; jm(n)=j; iuv(n)=2;
        N(i,j,2)=n; % stores indices
      end
  end
end

% only fill spots that dont correspond to a masked area with val, rest has value of 0
% odd indices are u vector
% even index v vector

% periodic boundary conditions:
N(1,:)=N(par.ni-1,:);
N(par.ni,:)=N(2,:);

neqns=n;

% initialize setup of sparse matrix:
nzmax=neqns*20;
nA=0;
i_A=zeros(nzmax,1);
j_A=zeros(nzmax,1);
v_A=zeros(nzmax,1);

RHS=zeros(neqns,1);

%% setup matrix and rhs: loop over grid points (rows of A), and set columns according to coeff of u and v at each grid point:
for i=2:par.ni-1
  for j=3:par.nj-2
    if mask(i,j)==1 % loop only through elements without mask
      if (mask(i+1,j)==1 ||  mask(i-1,j)==1) && (mask(i,j+1)==1 ||  mask(i,j-1)==1)
          B_iphalf_j=(B(i,j)+B(i+1,j))/2;
          B_imhalf_j=(B(i,j)+B(i-1,j))/2;
          B_i_jphalf=(B(i,j)+B(i,j+1))/2;
          B_i_jmhalf=(B(i,j)+B(i,j-1))/2;
    
          % sine at half locations:
          s_jphalf=(par.s(j)+par.s(j+1))/2;
          s_jmhalf=(par.s(j)+par.s(j-1))/2;
    
          %% u-equations:
          irow=N(i,j,1);
          n=N(i,j-1,1);
          if n~=0
              nA=nA+1; i_A(nA)=irow; j_A(nA)=n; 
              v_A(nA)=(0.5*B_i_jmhalf*(s_jmhalf^2/par.s(j-1))/par.dtheta_rad^2 ...
                   -mask(i,j-1)*mask(i,j+1)*(par.c(j)/par.s(j))*B(i,j)*0.5*par.s(j)^2/(par.s(j-1)*2*par.dtheta_rad));
          end
          n=N(i-1,j,1);
          if n~=0
              nA=nA+1; i_A(nA)=irow; j_A(nA)=n; 
              v_A(nA)=(2*B_imhalf_j*(1/par.s(j))/par.dphi_rad^2);
          end
          n=N(i,j,1);
          if n~=0
              nA=nA+1; i_A(nA)=irow; j_A(nA)=n; 
              v_A(nA)=(-2*B_iphalf_j*(1/par.s(j))/par.dphi_rad^2-2*B_imhalf_j* ...
                  (1/par.s(j))/par.dphi_rad^2-0.5*B_i_jphalf*(s_jphalf^2/par.s(j))/ ...
                  par.dtheta_rad^2-0.5*B_i_jmhalf*(s_jmhalf^2/par.s(j))/par.dtheta_rad^2);
          end
          n=N(i+1,j  ,1);
          if n~=0
              nA=nA+1;i_A(nA)=irow;j_A(nA)=n;
              v_A(nA)=(2*B_iphalf_j*(1/par.s(j))/par.dphi_rad^2);
          end
          n=N(i  ,j+1,1);
          if n~=0
              nA=nA+1;i_A(nA)=irow;j_A(nA)=n;
              v_A(nA)=(0.5*B_i_jphalf*(s_jphalf^2/par.s(j+1))/par.dtheta_rad^2+mask(i,j-1)* ...
                  mask(i,j+1)*(par.c(j)/par.s(j))*B(i,j)*0.5*par.s(j)^2/(par.s(j+1)*2*par.dtheta_rad));
          end
          n=N(i-1,j-1,2);
          if n~=0;nA=nA+1;i_A(nA)=irow;j_A(nA)=n;v_A(nA)=(0.5*B(i,j-1)/(4*par.dphi_rad*par.dtheta_rad) ...
                  +B(i-1,j)/(4*par.dphi_rad*par.dtheta_rad) );
          end
          n=N(i+1,j-1,2);
          if n~=0
              nA=nA+1;i_A(nA)=irow;j_A(nA)=n;
              v_A(nA)=(-0.5*B(i,j-1)/(4*par.dphi_rad*par.dtheta_rad)-B(i+1,j)/(4*par.dphi_rad*par.dtheta_rad) );
          end       
          n=N(i-1,j  ,2);
          if n~=0
              nA=nA+1;i_A(nA)=irow;j_A(nA)=n;v_A(nA)=(-2*(par.c(j)/par.s(j))*B(i-1,j)/(2*par.dphi_rad) ...
                  -mask(i-1,j)*mask(i+1,j)*(par.c(j)/par.s(j))*B(i,j)*0.5/(2*par.dphi_rad));
          end
          n=N(i+1,j  ,2);
          if n~=0
              nA=nA+1; i_A(nA)=irow; j_A(nA)=n; v_A(nA)=(2*(par.c(j)/par.s(j))*B(i+1,j)/(2*par.dphi_rad) ...
                  +mask(i-1,j)*mask(i+1,j)*(par.c(j)/par.s(j))*B(i,j)*0.5/(2*par.dphi_rad));
          end
          n=N(i-1,j+1,2);
          if n~=0
              nA=nA+1;
              i_A(nA)=irow;
              j_A(nA)=n;
              v_A(nA)=(-0.5*B(i,j+1)/(4*par.dphi_rad*par.dtheta_rad)-B(i-1,j)/(4*par.dphi_rad*par.dtheta_rad) );
          end
          n=N(i+1,j+1,2);
          if n~=0
              nA=nA+1; i_A(nA)=irow; j_A(nA)=n; v_A(nA)=(0.5*B(i,j+1)/(4*par.dphi_rad*par.dtheta_rad) ...
                  +B(i+1,j)/(4*par.dphi_rad*par.dtheta_rad) );
          end
          
          %% rhs for u equation:
          RHS(irow)= par.g*par.rho_I*h(i,j)*(1-par.mu)*( ...
              mask(i+1,j)*(h(i+1,j)-h(i,j))+mask(i-1,j)*(h(i,j)-h(i-1,j)) ...
              )/((mask(i+1,j)+mask(i-1,j))*par.dphi_rad);
    
          %% v-equations:
          irow=N(i,j,2);
          n=N(i-1,j-1,1);
          if n~=0
              nA=nA+1;i_A(nA)=irow;j_A(nA)=n;v_A(nA)=(0.5*B(i-1,j)*(par.s(j)/par.s(j-1))/(4*par.dphi_rad*par.dtheta_rad) ...
                  +par.s(j)*B(i,j-1)*(1/par.s(j-1))/(4*par.dphi_rad*par.dtheta_rad));
          end
          n=N(i+1,j-1,1);
          if n~=0
              nA=nA+1;i_A(nA)=irow;j_A(nA)=n;v_A(nA)=(-0.5*B(i+1,j)*(par.s(j)/par.s(j-1))/(4*par.dphi_rad*par.dtheta_rad) ...
                  -par.s(j)*B(i,j-1)*(1/par.s(j-1))/(4*par.dphi_rad*par.dtheta_rad) );
          end
          n=N(i-1,j  ,1);
          if n~=0
              nA=nA+1;i_A(nA)=irow;j_A(nA)=n;v_A(nA)=(+mask(i-1,j)*mask(i+1,j)*(par.c(j)/par.s(j))*B(i,j)/(2*par.dphi_rad));
          end
          n=N(i+1,j  ,1);
          if n~=0;nA=nA+1;i_A(nA)=irow;j_A(nA)=n;v_A(nA)=(-mask(i-1,j)*mask(i+1,j)*(par.c(j)/par.s(j))*B(i,j)/(2*par.dphi_rad));
          end
          n=N(i-1,j+1,1);
          if n~=0;nA=nA+1;i_A(nA)=irow;j_A(nA)=n;v_A(nA)=(-0.5*B(i-1,j)*(par.s(j)/par.s(j+1))/(4*par.dphi_rad*par.dtheta_rad) ...
                  -par.s(j)*B(i,j+1)*(1/par.s(j+1))/(4*par.dphi_rad*par.dtheta_rad));
          end
          n=N(i+1,j+1,1);
          if n~=0;nA=nA+1;i_A(nA)=irow;j_A(nA)=n;v_A(nA)=(0.5*B(i+1,j)*(par.s(j)/par.s(j+1))/(4*par.dphi_rad*par.dtheta_rad) ...
                                                        +par.s(j)*B(i,j+1)*(1/par.s(j+1))/(4*par.dphi_rad*par.dtheta_rad));
          end
          n=N(i,j-1,2);
          if n~=0
              nA=nA+1;i_A(nA)=irow;j_A(nA)=n;v_A(nA)=(B_i_jmhalf*s_jmhalf/par.dtheta_rad^2 ...
                                                        +par.s(j)*B_i_jmhalf*(1/s_jmhalf)*par.s(j-1)/par.dtheta_rad^2);
          end
          n=N(i-1,j,2);
          if n~=0
              nA=nA+1;i_A(nA)=irow;j_A(nA)=n;v_A(nA)=(0.5*B_imhalf_j*(1/par.s(j))/par.dphi_rad^2);
          end
          n=N(i,j,2);
          if n~=0
              nA=nA+1;i_A(nA)=irow;j_A(nA)=n;v_A(nA)=(-0.5*B_iphalf_j*(1/par.s(j))/par.dphi_rad^2 ...
                  -0.5*B_imhalf_j*(1/par.s(j))/par.dphi_rad^2-B_i_jphalf*s_jphalf/par.dtheta_rad^2 ...
                  -B_i_jmhalf*s_jmhalf/par.dtheta_rad^2-par.s(j)^2*B_i_jmhalf*(1/s_jmhalf)/par.dtheta_rad^2 ...
                  -par.s(j)^2*B_i_jphalf*(1/s_jphalf)/par.dtheta_rad^2-(par.c(j)/par.s(j))*B(i,j)*par.c(j));
          end
          n=N(i+1,j,2);
          if n~=0
              nA=nA+1;i_A(nA)=irow;j_A(nA)=n;v_A(nA)=(0.5*B_iphalf_j*(1/par.s(j))/par.dphi_rad^2);
          end
          n=N(i,j+1,2);
          if n~=0
              nA=nA+1;
              i_A(nA)=irow;
              j_A(nA)=n;
              v_A(nA)=(B_i_jphalf*s_jphalf/par.dtheta_rad^2 ...
                  +par.s(j)*B_i_jphalf*(1/s_jphalf)*par.s(j+1)/par.dtheta_rad^2);
          end
          
          %% rhs for v equation:
           RHS(irow)= par.s(j)*par.g*par.rho_I*(1-par.mu)*h(i,j)*(...
                 mask(i,j+1)*(h(i,j+1)-h(i,j))+mask(i,j-1)*(h(i,j)-h(i,j-1)) ...
                 )/((mask(i,j+1)+mask(i,j-1))*par.dtheta_rad);
           % RHS = 0 outside of if struct
      else
          RHS(irow) = 0;
      end
    end
  end
end

A=sparse(i_A(1:nA),j_A(1:nA),v_A(1:nA),neqns,neqns);

if nnz(A)>=nzmax
  % test consistency of sparse representation of A:
  fprintf(1,' *** too many nonzero elements in A: nnz(A)=%d, nzmax=%d\n',nnz(A),nzmax);
  pause
end

if 0
  figure(54); clf
  fprintf(1,'function solve_uv_2d_matrix_form: plotting A, N & RHS...');
  subplot(2,2,1); imagesc(A(1:neqns,1:neqns)/max(max(A))); title('A'); colorbar;
  subplot(2,2,2); imagesc(squeeze(N(:,:,1)'))/max(max(max(N(:,:,1)))); title('N(:,:,1)'); colorbar;
  subplot(2,2,3); imagesc(squeeze(N(:,:,2)'))/max(max(max(N(:,:,2)))); title('N(:,:,2)'); colorbar;
  subplot(2,2,4); plot(RHS); title('RHS');
  fprintf(1,' done.\n');
  pause
end


% x = A\B is a solution to the equation A*x = B, if it exists. for x=UV
ind = ~any(A~=0, 2);
A(ind, ind) = 0.5*speye(nnz(ind));
UV=A\RHS;
check_nan(RHS); % check for inconsistencies

%% put back on i,j grid:
u=zeros(par.ni,par.nj); v=u; 
plot_RHS=0; % use to debug
if plot_RHS; RHS_u=v; RHS_v=v; end

for n=1:neqns
  if iuv(n)==1
    u(im(n),jm(n))=UV(n);
    if plot_RHS 
        RHS_u(im(n),jm(n))=RHS(n); 
    end
  elseif iuv(n)==2
    v(im(n),jm(n))=UV(n);
    if plot_RHS 
        RHS_v(im(n),jm(n))=RHS(n); 
    end
  end
end

if plot_RHS
  figure(55); clf; 
  subplot(1,2,1); imagesc(RHS_u'); title(sprintf('RHS_u, iter_eff_viscosity=%g',par.iter_eff_viscosity)); colorbar; 
  subplot(1,2,2); imagesc(RHS_v'); title('RHS_v'); colorbar; 
  pause
end


%% north and south b.c.:
u(:,2)=u(:,2)*0;
u(:,par.nj-1)=v(:,par.nj-1)*0;
v(:,2)=v(:,2)*0;
v(:,par.nj-1)=v(:,par.nj-1)*0;
u(:,1)=u(:,1)*0;
u(:,par.nj)=v(:,par.nj)*0;
v(:,1)=v(:,1)*0;
v(:,par.nj)=v(:,par.nj)*0;

%% periodic boundary conditions:
u(1,:)=u(par.ni-1,:);
u(par.ni,:)=u(2,:);
v(1,:)=v(par.ni-1,:);
v(par.ni,:)=v(2,:);

if 0
  figure(55); clf; 
  contourf(u);
  colorbar;

  figure(56); clf; 
  contourf(v);
  colorbar;
  % plot terms in v-equation, try this when problem is indep of x:
  i=15;
  for j=2:par.nj-1
    B_i_jphalf=(B(i,j)+B(i,j+1))/2;
    B_i_jmhalf=(B(i,j)+B(i,j-1))/2;
    term_B(j)=(B_i_jphalf*(v(i,j+1)-v(i,j))-B_i_jmhalf*(v(i,j)-v(i,j-1)))/par.dtheta_rad^2;
    term_h(j)=par.g*par.rho_I*(1-par.mu)*(h(i,j+1)-h(i,j-1))/(2*par.dtheta_rad);
  end
  plot(term_B,'r','linewidth',2); hold on; plot(term_h,'--b'); title('terms in v eqn');
  legend('B','h');
  pause
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S=S_total(par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return upper-surface source-sink for ice as function of latitude.
% artificially divide pollard forcing as half top and half bottom:
load(par.S_filename);
if length(S_Pollard)~=par.nj
  disp('*** wrong dimensions of imported Pollard forcing.')
  fprintf(1,' filename: %s;\n par.nj=%d, length(S_Pollard)=%d\n',par.S_filename,par.nj,length(S_Pollard));
  pause
else
  S=S_Pollard;
end
% Convert from mm/year to use m/s:
S=S/(1000*par.year);
% ice accumulation -> positive, ice melt at lower latitude -> negative
% value is calculated by pollard from a 1D climate model, we don't use this
% value so for now set to 0, since it doesnt correspond to our GCM


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S=S_bot(par,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return basal melt/ freeze ice as function of latitude.
% sink -> perte en H due à la fonte

% This returns a 1D matrix of surface T, so needs to be modified (this
% takes the avg of the longitude (not applicable for tidal locking) and L
% is spatial resolution (this doesnt run, need to define the vars), for
% now we set this to zero, uncomment the last lines to use equations to 
% calc the bottom freeze / melt rate
 
S=0*h;
   
% T_s=T_surface_function(par);
% T_s = mean(T_s,2)'; adjust to 2D instead of 1d
% S=(-par.kappa_ice*(T_s-par.T_f)-par.F_g*h)./(par.rho_I*par.L*h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=A_g(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature dependence of ice viscosity, Goodman & Pierhumbert 2003.
if T<263.15
  A0=3.61e-13;
  Q=60e3;
else
  A0=1.734e3;
  Q=139e3;
end

R=8.3144621; % http://en.wikipedia.org/wiki/Gas_constant
A=A0*exp(-Q/(R*T));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TS=T_surface_function(par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y_middle=par.y(par.nj-1)/2;
lat=(par.y-y_middle)*0.5*pi/y_middle;

% if needed plot surface temperature profiles
if 0 && par.n == 1
    % longitunal and lateral profile for 'mine' and comparison
    fig_lat_temp = figure(10);clf;
    lat_arr = linspace(-90,90, 176);
    plot(lat_arr, par.tLat, '-b');
    hold on
    plot(lat_arr, par.T_f-79+48*cos(lat).^2, '-r');
    hold off
    xlim([-90,90]);
    title('Temperature Profiles');
    xlabel('Latitude (°)');
    ylabel('Temperature (k)');
    legend('Output','cold');
    savefig_pdf(par, fig_lat_temp, 'latitudinal_temp_profile');
    
    fig_lon_temp = figure(11);clf;
    lon_arr = linspace(0,360, 176);
    plot(lon_arr, par.tLon, '-b');
    xlim([0,360]);
    title('Temperature Profile');
    xlabel('Longitude (°)');
    ylabel('Temperature (k)'); 
    savefig_pdf(par, fig_lon_temp, 'longitudinal_temp_profile');
end 

if strcmp(par.T_surface_profile_type,'mine')
  TS = par.tSurf;
elseif strcmp(par.T_surface_profile_type,'cold')
  % Dorian's low CO2:
  temp=par.T_f-79+48*cos(lat).^2;
  TS = repmat(temp, 176, 1);

  fig_tsurf_cold = figure(12);clf;
  [~,hc] = contourf(par.phi,90-par.theta(end:-1:1),rot90(TS));
  set(hc,'EdgeColor','none');
  colormap('jet');
  colorbar;
  title('Surface Temperature (K)');
  xlabel('longitude','FontSize',20);
  ylabel('latitude','FontSize',20);
  savefig_pdf(par, fig_tsurf_cold, 'tsurf-cold');
elseif strcmp(par.T_surface_profile_type,'warm')
  % Dorian's high CO2:
  TS=par.T_f-58+47*cos(lat).^2;
else
  disp('*** no such surface temperature profile.');
  quit
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par = adjust_mask(par, h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modifies the dynamic mask according to the ice elevation that changes
par.ocean_mask(h<10) = 0; 
par.ocean_mask(h>10) = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_nan(mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for inconsistencies
if anynan(mat)
    if isvector(mat)
        row=find(isnan(mat));
        disp('problem!! NaN found');
        disp(row);
    else
        [row,col]=find(isnan(mat));
        disp('problem!! NaN found');
        disp([row,col]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savefig_pdf(par, handle, name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handle, 'PaperUnits', 'inches'); 
set(handle, 'PaperSize', [6 4]);
set(handle, 'PaperPosition', [0 0 6 4]); % [left, bottom, width, height];
saveas(handle,sprintf('Figures/fig-exp-%2.2d-%s.pdf',par.expnum, name));
%saveas(handle,sprintf('Figures/fig-exp-%2.2d-%s.fig',par.expnum, name));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B=calc_eff_viscosity(par,T_surface,h,u,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% effective viscosity for 2d model in spherical coordinates

if par.which_mask == 0
    mask=par.land_mask;
elseif par.which_mask == 1
    par = adjust_mask(par,h);
    mask=par.ocean_mask;
elseif par.which_mask == 2
    par = adjust_mask(par,h);
    mask=par.ocean_mask.*par.land_mask;
else
    mask=par.no_mask;
end

s=par.s; 
B=zeros(par.ni,par.nj);
dot_eps0=1e-14; 

% matrix 2D for T_surface and dependent linear temperature profile
% go line by line to fill out matrix of B
for i=2:par.ni-1
  for j=2:par.nj-1 

    % the eqn: [\dot\epsilon=sqrt(\dot\epsilon_ij\dot\epsilon_ij/2)]:
    if mask(i,j)== 0
        B(i,j)= 0; % only calculate dot_eps and B if there is no mask
    elseif (mask(i+1,j)==1 ||  mask(i-1,j)==1) && (mask(i,j+1)==1 ||  mask(i,j-1)==1)
        % create linear vertical temp profile for each val of T_surface
        
        % if surface temp has reached 273k -> vertical profile needs to be
        % done differently, for now arbitrary surface temp is set
        % if too warm -> set B(i,j)=0 and change mask
        if T_surface(i, j) > 272
            T=265:0.1:par.T_f;
        else
            T=T_surface(i, j):0.5:par.T_f; % for more precise approximation decrease incrementation (0.1)
        end
        
        % calculate temperature dependence of viscosity using Glen's flow
        % law
        AA=T*NaN; % create AA array same dim as T

        for k=1:length(T) 
            AA(k)=(A_g(T(k)))^(-1/3); % glens flow param
        end

        % compare the error from AA for testing approximations of glens
        % flow law
        if 0
            relDif = (mean(AA) - mean(AA_test))/mean(AA);
            if abs(relDif) > 0.05 && abs(relDif) ~= Inf % check if error is within 5%
                disp('problem');
                disp('stopping')
            end
        end 

        % taylor series for e^(-1/x) ????
        % in range 0 to 300
        % sum from 0 to inf [ (-1^n)/(n!) * x^-n ]
        % A=A0*exp(-Q/(R*T));
        % 1/e + (x-1)/e - (x-1)^2/2e - (x-1)^3/6e
        % A0=3.61e-13;
        % Q=60e3;
        % R=8.3144621;
        % factor = Q/R
        % AA_tt = 1/exp(factor) + 2(x-T) 

        % calculate dot_epsilon matrix elements
        eps_xx=(1/(par.R*s(j)))*( ...
                 (mask(i+1,j)*(u(i+1,j)-u(i,j))+mask(i-1,j)*(u(i,j)-u(i-1,j))) ...
                  /((mask(i+1,j)+mask(i-1,j))*par.dphi_rad) ...
                  +v(i,j)*par.c(j));
    
        eps_yy=(1/par.R)*...
                 (mask(i,j+1)*(v(i,j+1)-v(i,j))+mask(i,j-1)*(v(i,j)-v(i,j-1))) ...
                 /((mask(i,j+1)+mask(i,j-1))*par.dtheta_rad);
    
        eps_zz=-(eps_xx+eps_yy);
    
        eps_xy=(1/par.R)*0.5*( ...
            (1/s(j))*(mask(i+1,j)*(v(i+1,j)-v(i,j))+mask(i-1,j)*(v(i,j)-v(i-1,j)))...
            /((mask(i+1,j)+mask(i-1,j))*par.dphi_rad) ...
            +s(j)*(mask(i,j+1)*(u(i,j+1)/s(j+1)-u(i,j)/s(j))+mask(i,j-1)*(u(i,j)/s(j)-u(i,j-1)/s(j-1))) ...
            /((mask(i,j+1)+mask(i,j-1))*par.dtheta_rad));
    
         % eps_zx = 0.5*u(i,j)/h(i,j);
         % eps_zy = 0.5*v(i,j)/h(i,j);
         
         % eps_zx = du/dh % assume dz en x,y = 0
         % extra eps dot to account for grounding ? not sufficient
         % dot_eps=sqrt((eps_xx^2+eps_yy^2+eps_zz^2+2*eps_xy^2+2*eps_zx^2+2*eps_zy^2)/2);
         dot_eps=sqrt((eps_xx^2+eps_yy^2+eps_zz^2+2*eps_xy^2)/2);
 
         % calculate B with all the parameters
         B(i,j)=(h(i,j)/par.R)*mean(AA)*(dot_eps+dot_eps0)^(1/par.nn-1);
   
    else 
        % for values that have masked neighboring indices, we can't
        % calculate B, so set average value (could be a min instead)
        B(i,j)= par.B_ocean;
    end
  end
end

% check for inconsistencies
check_nan(B);
B(isnan(B))=par.B_ocean;

% periodic boundary conditions for B:
% adjusts last latitude to be equal to first latitude
B(1,:)=B(par.ni-1,:);
B(par.ni,:)=B(2,:);

