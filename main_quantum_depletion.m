% Quantum depletion
% DKS 06/02/2017

clear all; close all; clc;

%%% Constants
tof_vz=9.8*0.430;    % atom free-fall vert v at detector hit for T-to-Z conversion;

%%% USER INPUTS
path_config='C:\Users\HE BEC\Documents\MATLAB\quantum-depletion\config_060217_test.m';

% vars to save to output
vars_save={'path_config',...
    'zxy_0','files_out',...
    'k_ff'};


%% Main
% load config
run(path_config);

% Load TXY data and crop to region of interest
[txy_raw,files_out]=loadExpData(configs,verbose);
nShot=size(txy_raw,1);

%%% Pre-processing raw data
txy_0=cell(nShot,1);    % oscillation compensated txy counts
zxy_0=cell(nShot,1);    % T-Z conversion
for i=1:nShot
    cond_cent_tmp=mean(txy_raw{i},1);   % approx condensate centre from average of captured
    n_count_tmp=size(txy_raw{i},1);     % number of counts in this shot
    txy_0{i}=txy_raw{i}-repmat(cond_cent_tmp,[n_count_tmp,1]);  % centre around self-average
    
    % T-Z conversion
    zxy_0{i}=txy_0{i};
    zxy_0{i}(:,1)=zxy_0{i}(:,1)*tof_vz;    % TOF - dz at time of detection
end

% Plot far-field ZXY
if verbose>1    
    h_zxy_ff=figure();
    plot_zxy(zxy_0,100,'r');
    title('Condensate point cloud');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    view(3);
    axis equal;
    
    % save plot
    fname_str='zxy_ff';
    saveas(h_zxy_ff,[configs.files.dirout,fname_str,'.png']);
    saveas(h_zxy_ff,[configs.files.dirout,fname_str,'.fig']);
end


%%% far-field momentum distribution
k_ff=abs(vertcat(zxy_0{:}));     % collate all shots and convert to abs(k)-space (TODO)

% Radial (ZY plane)
k_rad=k_ff(:,[1,3]);
log_k_rad=real(log(k_rad));   % get log dist

[n_hist_krad,ed_hist_krad]=histcounts(k_rad,100);       % lin hist
[n_hist_lgkrad,ed_hist_lgkrad]=histcounts(log_k_rad,100);   % log hist

% Longitudinal (X-axis)
k_lon=k_ff(:,2);
log_k_lon=real(log(k_lon));   % get log dist

[n_hist_klon,ed_hist_klon]=histcounts(k_lon,100);   % lin hist
[n_hist_lgklon,ed_hist_lgklon]=histcounts(log_k_lon,100);  % log hist

% Plot
%%% Linear-k-ff histogram
h_hist_k=figure();
hold on;
plot((ed_hist_krad(1:end-1)+ed_hist_krad(2:end))/2,...
    n_hist_krad,'*--');
plot((ed_hist_klon(1:end-1)+ed_hist_klon(2:end))/2,...
    n_hist_klon,'o-');
legend({'Radial','Longitudinal'});

%%% Llog-k-ff histogram
h_hist_lgk=figure();
loglog((ed_hist_lgkrad(1:end-1)+ed_hist_lgkrad(2:end))/2,...
    n_hist_lgkrad,'*--');
hold on;
loglog((ed_hist_lgklon(1:end-1)+ed_hist_lgklon(2:end))/2,...
    n_hist_lgklon,'o-');
legend({'Radial','Longitudinal'});


%% Save data
for i = 1:length(vars_save)
    if ~exist(vars_save{i},'var')
        warning(['Variable "',vars_save{i},'" does not exist.']);
        continue;
    end
    save(configs.files.saveddata,vars_save{i},'-append');     % to overcome version conflict
end