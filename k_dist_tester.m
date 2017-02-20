%%%% Tests k-distribution algorithm against 3D flat background counts
vz=0.416*9.8;   % vz at tof of detection


run('flat_background.m');   % builds TXY_FLAT - a noisy experimental data simulator

%% T to Z
zxy_flat=TXY_FLAT;
nShot=size(TXY_FLAT,1);
for iShot=1:nShot
    zxy_flat{iShot}(:,1)=vz*zxy_flat{iShot}(:,1);
end

% display pointcloud
n_shot_plot=20;
if n_shot_plot>nShot
    n_shot_plot=nShot;
end
h_pointcloud=figure();
plot_zxy(zxy_flat(1:n_shot_plot));
view(3);
axis equal;
xlabel('X');ylabel('Y');zlabel('Z');

zxy_all=vertcat(zxy_flat{:});

%% Cartesian to cylindrical coordinate conversion
% [R_YZ, theta, Perpendicular height]
r_cyl_flat=cell(nShot,1);
for iShot=1:nShot
    % get cartesian coords
    Z=zxy_flat{iShot}(:,1);
    X=zxy_flat{iShot}(:,2);
    Y=zxy_flat{iShot}(:,3);
    
    % calculate cylindrical coords
    R_yz=sqrt(Z.^2+Y.^2);
    H_perp=X;
    Theta=atan2(Z,Y)+pi/2;  % define -Z to be 0 and measure by RH +X
    
    % store cyl coords
    r_cyl_flat{iShot}=[R_yz,Theta,H_perp];
end

% collate all shots
r_cyl_all=vertcat(r_cyl_flat{:});

%% Cylindrical slice
r_cyl_slice=r_cyl_all;

% define cyl-slice region
cyl_theta_lims=pi+deg2rad(15)*[-0.5,0.5];  % azimuthal angle
cyl_rad_lims=[0,0.3];    % radial disp from origin
cyl_trans_lims=0.02*[-0.5,0.5];  % transverse displacement from origin

% cull in each dimension
idx_all=1:size(r_cyl_all,1);

% 1. angular
theta_temp=r_cyl_all(:,2);
theta_temp=wrapTo2Pi(theta_temp-cyl_theta_lims(1));     % zero around theta_i
idx_in=(theta_temp<diff(cyl_theta_lims));   % get indices which are in angular cropping region
idx_all=idx_all(idx_in);        % update all index
r_cyl_slice=r_cyl_slice(idx_in,:);    % cull in theta

% 2. radial
rad_temp=r_cyl_slice(:,1);
idx_in=(rad_temp>cyl_rad_lims(1))&(rad_temp<cyl_rad_lims(2));
idx_all=idx_all(idx_in);
r_cyl_slice=r_cyl_slice(idx_in,:);

% 3. transverse
trans_temp=r_cyl_slice(:,3);
idx_in=(trans_temp>cyl_trans_lims(1))&(trans_temp<cyl_trans_lims(2));
idx_all=idx_all(idx_in);
r_cyl_slice=r_cyl_slice(idx_in,:);

% Display results
h_slice_capture=figure();
hold on;
scatter_zxy(zxy_all);
scatter_zxy(zxy_all(idx_all,:),30,'r');
view(3);
axis equal;
xlabel('X');ylabel('Y');zlabel('Z');

%% Radial density profile in slice
% get 1D data
r_cyl=r_cyl_slice(:,1);

% Log-spaced bin properties
nbin_log=100;
r_log_edge=logspace(log10(0.01),log10(0.3),nbin_log);             % log spaced edges
r_log_cent=sqrt(r_log_edge(1:end-1).*r_log_edge(2:end));    % geometric avgd centre

% bin volume
bin_volume=diff(cyl_trans_lims)*diff(cyl_theta_lims)*r_log_cent.*diff(r_log_edge);

% histogram
N_log=histcounts(r_cyl,r_log_edge);

% number density
nden_log=N_log./(nShot*bin_volume); % no QE - very basic

% test flat density
vol=vz*diff(LIM{1})*diff(LIM{2})*diff(LIM{3});     % volume [m/^3]
nden_flat=N_PER_SHOT/vol;   % average flat count density [/m^3]

% Plot resulting histogram and density profile
h_r_profile=figure();
subplot(1,2,1);
barwidth=1;
plot(r_log_cent,N_log,'.');

subplot(1,2,2);
plot(r_log_cent,nden_log,'*');
set(gca,'xScale','log');
set(gca,'yScale','log');