% Quantum depletion
% DKS 06/02/2017

clear all; close all; clc;

%%% USER INPUTS
path_config='C:\Users\HE BEC\Documents\MATLAB\quantum-depletion\config_060217_test.m';

% vars to save to output
vars_save={'path_config',...
    'zxy_0','files_out',...
    'k_ff',...
    'N_hist_krad','ed_hist_krad',...
    'N_hist_lgkrad','ed_hist_lgkrad',...
    'N_hist_klon','ed_hist_klon',...
    'N_hist_lgklon','ed_hist_lgklon',...
    'n_hist_krad','n_hist_klon',...
    'zxy_slice','r_1D',...
    'N_r1D','ed_r1D',...
    'area_perp_1D',...
    'n_r1D','n_kff',...
    'log_r_1D','N_lgr1D','ed_lgr1D',...
    'dr1D','n_lgr1D','n_kff_lg',...
    'ed_kff','ed_kff_log',...
    };


%% Main
% load config
run(path_config);

% Load misc params
hbar=configs.const.hbar;
m_He=configs.const.m_He;
tof=configs.const.tof;
vz=9.8*tof;     % atom free-fall vert v at detector hit for T-to-Z conversion;
detQE=configs.const.detect_qe;

% Load TXY data and crop to region of interest
[txy_raw,files_out]=loadExpData(configs,verbose);
nShot_raw=size(txy_raw,1);


%% Pre-processing
% filter shots with spurious counts in 1D window of interest
txy_0=cell(nShot_raw,1);        % oscillation compensated txy counts
zxy_0=cell(nShot_raw,1);        % T-Z conversion
zxy_slice=cell(nShot_raw,1);    % 1D slice captured counts

Nmin_1D=configs.slice.mincount;     % minimum count in 1D slice to pass
N_1D_all=zeros(nShot_raw,1);        % number captured in 1D slice

fid_analysis=zeros(nShot_raw,1);       % fild id's used for following analysis (post-filter)
counter=1;
for i=1:nShot_raw
    cond_cent_tmp=mean(txy_raw{i},1);   % approx condensate centre from average of captured
    n_count_tmp=size(txy_raw{i},1);     % number of counts in this shot
    txy_0_temp=txy_raw{i}-repmat(cond_cent_tmp,[n_count_tmp,1]);  % centre around self-average
    
    % T-Z conversion
    zxy_0_temp=txy_0_temp;
    zxy_0_temp(:,1)=zxy_0_temp(:,1)*vz;    % TOF - dz at time of detection

    % Take 1D slice
    [zxy_slice_temp,~,n_captured]=cylindercull(zxy_0_temp,configs.slice.cyl_cent,...
        configs.slice.cyl_dim,configs.slice.cyl_orient);
    N_1D_all(i)=n_captured;         % save captured number
    
    % filter shots based on counts in 1D slice
    if n_captured>Nmin_1D
        txy_0{counter}=txy_0_temp;
        zxy_0{counter}=zxy_0_temp;
        zxy_slice{counter}=zxy_slice_temp;
        fid_analysis(counter)=files_out.id_ok(i);
        
        counter=counter+1;
    else
        if verbose>0
            warning('1D filter: Low-count in file #%d. Discarding from further processing.',files_out.id_ok(i));
        end
    end
end
nShot=counter-1;    % number of shots passed from filtering for N in 1D slice
% clean up arrays
txy_0=txy_0(1:nShot);
zxy_0=zxy_0(1:nShot);
zxy_slice=zxy_slice(1:nShot);
fid_analysis=fid_analysis(1:nShot);

%%% PLOT
% Plot far-field ZXY
if verbose>2    
    h_zxy_ff=figure();
    plot_zxy(zxy_0,20,'k');
    title('Condensate point cloud');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    view(3);
    axis equal;
    
    % save plot
    fname_str='zxy_ff';
    saveas(h_zxy_ff,[configs.files.dirout,fname_str,'.png']);
    saveas(h_zxy_ff,[configs.files.dirout,fname_str,'.fig']);
end

% Plot 1D slice point cloud
if verbose>1
    h_zxy_slice=figure();
    plot_zxy(zxy_slice,100,'k');
    title('1D sampled condensate point cloud');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    view(3);
    axis equal;
    
    % save plot
    fname_str='zxy_slice';
    saveas(h_zxy_slice,[configs.files.dirout,fname_str,'.png']);
    saveas(h_zxy_slice,[configs.files.dirout,fname_str,'.fig']);
end

%% Plot summary
h_zxy_summ=figure();
nShotSumm=50;   % number of shots to plot as summary
if nShot<nShotSumm
    nShotSumm=nShot;
end
plot_zxy(zxy_0(1:nShotSumm),10,'r');
hold on;
plot_zxy(zxy_slice(1:nShotSumm),300,'b');

title('Summary - Condensate point cloud');
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
view(3);
axis equal;

% save plot
fname_str='zxy_summ';
saveas(h_zxy_summ,[configs.files.dirout,fname_str,'.png']);
saveas(h_zxy_summ,[configs.files.dirout,fname_str,'.fig']);


%% Density profiling
% collate and collapse all captured count to 1D
r_1D=abs(vertcat(zxy_slice{:}));
r_1D=r_1D(:,configs.slice.cyl_orient);      % 1D squeezed data in real-space

% r to k
k_1D=r2k(r_1D);     % array of 1D k in [m^-1]

%%% LINEAR DIST
hist_nbin=configs.hist.nbin;
[N_r1D,ed_r1D]=histcounts(r_1D,hist_nbin);  % lin hist - autoscale bins to lim
N_r1D=N_r1D/(nShot*detQE);  % normalise for a single shot and consider detector QE

% number to density
area_perp_1D=pi*(configs.slice.cyl_rad^2);  % area of integrated dims
n_r1D=N_r1D./(diff(ed_r1D)*area_perp_1D);   % number density [m^-3]

n_kff=(hbar*tof/m_He)^3*n_r1D;  % number density (real space) to far-field momentum density


if verbose>0	% plot
    % real space density dist
    h_nr1D=figure();
    plot(1e3*(ed_r1D(1:end-1)+ed_r1D(2:end))/2,...
        1e-9*n_r1D,'--');
    title('1D condensate number profile');
    xlabel('$r$ [mm]'); ylabel('$n(r,\overline{t})$ [mm$^{-3}$]');
    
    fname_str='nr1D';
    saveas(h_nr1D,[configs.files.dirout,fname_str,'.png']);
    saveas(h_nr1D,[configs.files.dirout,fname_str,'.fig']);
    
    % far-field momentum space (lin)
    h_nkff=figure();
    ed_kff=(m_He/(hbar*tof))*ed_r1D;
    plot(1e-6*(ed_kff(1:end-1)+ed_kff(2:end))/2,...
        1e18*n_kff,'-');
    title('1D condensate momentum profile');
    xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
    
    fname_str='nkff';
    saveas(h_nkff,[configs.files.dirout,fname_str,'.png']);
    saveas(h_nkff,[configs.files.dirout,fname_str,'.fig']);
end

%%% LOG DIST (Method 1)
log_r_1D=real(log(r_1D));   % convert to log space
[N_lgr1D,ed_lgr1D]=histcounts(log_r_1D,hist_nbin);
N_lgr1D=N_lgr1D/(nShot*detQE);      % normalise for single shot and detector QE
dr1D=exp(ed_lgr1D(2:end))-exp(ed_lgr1D(1:end-1));   % hist bin width in 1D
n_lgr1D=N_lgr1D./(dr1D*area_perp_1D);       % number density [m^-3]

n_kff_lg=(hbar*tof/m_He)^3*n_lgr1D; % ff momentum density

if verbose>0    % plot
    % far-field momentum space (log)
    h_nkff_log=figure();
    ed_kff_log=(m_He/(hbar*tof))*exp(ed_lgr1D);
    loglog(1e-6*(ed_kff_log(1:end-1)+ed_kff_log(2:end))/2,...
        1e18*n_kff_lg,'*-');
    xlim([0.1,10]);   %   limit x-axis to like Clement paper
    grid on;
    title('1D condensate momentum profile');
    xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
    
    fname_str='nkff_log';
    saveas(h_nkff_log,[configs.files.dirout,fname_str,'.png']);
    saveas(h_nkff_log,[configs.files.dirout,fname_str,'.fig']);
end

%%% LOG DIST (Method 2)
% TODO specify bin edges for log-space
% ed_lgk=configs.hist.ed_lgk;     % get log-spaced edges
% N_lgr1D_test=histcounts(k_1D,ed_lgk);

%%% n(k)k4 scaled plot
kcent=(ed_kff_log(1:end-1)+ed_kff_log(2:end))/2;
n_scaled=n_kff_lg.*(kcent.^4);

if verbose>0    % plot
    h_scaled_nk=figure();
    semilogy(1e-6*kcent,n_scaled,'*');
    ylim([1e8,1e11]);
    xlabel('$k$ [$\mu$m$^{-1}$]');
    ylabel('$k^{4}n_{\infty}(k)$ [m$^{-1}$]');
    
    fname_str='k4_scaled_dist';
    saveas(h_scaled_nk,[configs.files.dirout,fname_str,'.png']);
    saveas(h_scaled_nk,[configs.files.dirout,fname_str,'.fig']);
end

%% Fit to density profile



%% far-field momentum distribution
% k_ff=abs(vertcat(zxy_0{:}));     % collate all shots and convert to abs(k)-space (TODO)
% 
% % Radial (ZY plane)
% k_rad=k_ff(:,[1,3]);
% log_k_rad=real(log(k_rad));   % get log dist
% 
% [N_hist_krad,ed_hist_krad]=histcounts(k_rad,hist_nbin);       % lin hist
% [N_hist_lgkrad,ed_hist_lgkrad]=histcounts(log_k_rad,hist_nbin);   % log hist
% 
% % Longitudinal (X-axis)
% k_lon=k_ff(:,2);
% log_k_lon=real(log(k_lon));   % get log dist
% 
% [N_hist_klon,ed_hist_klon]=histcounts(k_lon,hist_nbin);   % lin hist
% [N_hist_lgklon,ed_hist_lgklon]=histcounts(log_k_lon,hist_nbin);  % log hist
% 
% %%% number to density conversion
% % TODO: scale to 1D
% % linear bin
% n_hist_krad=N_hist_krad./diff(ed_hist_krad);    %currently [m-1]
% n_hist_klon=N_hist_klon./diff(ed_hist_klon);
% 
% % logarithmic bin
% n_hist_lgkrad=N_hist_lgkrad./(exp(ed_hist_lgkrad(2:end))-exp(ed_hist_lgkrad(1:end-1)));
% n_hist_lgklon=N_hist_lgklon./(exp(ed_hist_lgklon(2:end))-exp(ed_hist_lgklon(1:end-1)));
% 
% % Plot
% %%% Linear r histogram
% h_hist_r=figure();
% hold on;
% plot((ed_hist_krad(1:end-1)+ed_hist_krad(2:end))/2,...
%     N_hist_krad,'--');
% plot((ed_hist_klon(1:end-1)+ed_hist_klon(2:end))/2,...
%     N_hist_klon,'-');
% 
% legend({'Radial','Longitudinal'});
% title('1D condensate number profile (linear)');
% xlabel('$r$ [m]'); ylabel('$N(r)$');
% 
% % save plot
% fname_str='hist_k';
% saveas(h_hist_r,[configs.files.dirout,fname_str,'.png']);
% saveas(h_hist_r,[configs.files.dirout,fname_str,'.fig']);
% 
% %%% linear r density
% h_nr=figure();
% hold on;
% plot((ed_hist_krad(1:end-1)+ed_hist_krad(2:end))/2,...
%     n_hist_krad,'--');
% plot((ed_hist_klon(1:end-1)+ed_hist_klon(2:end))/2,...
%     n_hist_klon,'-');
% 
% legend({'Radial','Longitudinal'});
% title('1D condensate density profile (linear)');
% xlabel('$r$ [m$^{-1}$]'); ylabel('$n(r)$');
% 
% % save plot
% fname_str='nr';
% saveas(h_nr,[configs.files.dirout,fname_str,'.png']);
% saveas(h_nr,[configs.files.dirout,fname_str,'.fig']);
% 
% 
% %%% loglog r histogram
% % Radial
% h_hist_lgkrad=figure();
% loglog((ed_hist_lgkrad(1:end-1)+ed_hist_lgkrad(2:end))/2,...
%     N_hist_lgkrad,'--');
% 
% title('1D radial condensate number profile (log-log)');
% xlabel('$r$ [m]'); ylabel('$N(r)$');
% axis tight;
% 
% % save plot
% fname_str='hist_lgkrad';
% saveas(h_hist_lgkrad,[configs.files.dirout,fname_str,'.png']);
% saveas(h_hist_lgkrad,[configs.files.dirout,fname_str,'.fig']);
% 
% % Longitudinal
% h_hist_lgklon=figure();
% loglog((ed_hist_lgklon(1:end-1)+ed_hist_lgklon(2:end))/2,...
%     N_hist_lgklon,'-');
% 
% title('1D longitudinal condensate number profile (log-log)');
% xlabel('$r$ [m]'); ylabel('$N(r)$');
% axis tight;
% 
% % save plot
% fname_str='hist_lgklon';
% saveas(h_hist_lgklon,[configs.files.dirout,fname_str,'.png']);
% saveas(h_hist_lgklon,[configs.files.dirout,fname_str,'.fig']);
% 
% %%% log r density
% % Radial
% h_n_lgkrad=figure();
% loglog((ed_hist_lgkrad(1:end-1)+ed_hist_lgkrad(2:end))/2,...
%     n_hist_lgkrad,'--');
% 
% title('1D radial condensate density profile');
% xlabel('$r$ [m]'); ylabel('$n(r)$ [m$^{-1}$]');
% axis tight;
% 
% % save plot
% fname_str='n_lgkrad';
% saveas(h_n_lgkrad,[configs.files.dirout,fname_str,'.png']);
% saveas(h_n_lgkrad,[configs.files.dirout,fname_str,'.fig']);
% 
% % Longitudinal
% h_n_lgklon=figure();
% loglog((ed_hist_lgklon(1:end-1)+ed_hist_lgklon(2:end))/2,...
%     n_hist_lgklon,'-');
% 
% title('1D longitudinal condensate density profile');
% xlabel('$r$ [m]'); ylabel('$n(r)$ [m$^{-1}$]');
% axis tight;
% 
% % save plot
% fname_str='n_lgklon';
% saveas(h_n_lgklon,[configs.files.dirout,fname_str,'.png']);
% saveas(h_n_lgklon,[configs.files.dirout,fname_str,'.fig']);


%% Save data
for i = 1:length(vars_save)
    if ~exist(vars_save{i},'var')
        warning(['Variable "',vars_save{i},'" does not exist.']);
        continue;
    end
    save(configs.files.saveddata,vars_save{i},'-append');     % to overcome version conflict
end