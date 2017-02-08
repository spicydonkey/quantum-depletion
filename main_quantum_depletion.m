% Quantum depletion
% DKS 06/02/2017

clear all; close all; clc;

%%% USER INPUTS
path_config='C:\Users\HE BEC\Documents\MATLAB\quantum-depletion\config_060217_test.m';

% vars to save to output
vars_save={'path_config',...
    'zxy_0','files_out',...
    'zxy_slice',...
    'r_1D','k_1D',...
    'r_perp_area','k_perp_area'...
    'hist_r1D','nden_r1D','nden_k1D',...
    'hist_lgk1D','nden_lgk1D',...
    'nk4_scaled',...
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

%% Proprocessing summary
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


%% 1D cull and processing
% collate and collapse all captured count to 1D
r_1D=abs(vertcat(zxy_slice{:}));
r_1D=r_1D(:,configs.slice.cyl_orient);	% 1D culled data in real-space

% r to k
k_1D=r2k(r_1D);     % array of 1D k in [m^-1]


%% Density profiling
r_perp_area=pi*(configs.slice.cyl_rad^2);       % real area of integrated dims
k_perp_area=pi*(r2k(configs.slice.cyl_rad)^2);  % k area int area

%%% LINEAR DIST
hist_r1D.binN=configs.hist.nbin;    % get number of bins to set up auto
[hist_r1D.N,hist_r1D.binEdge]=histcounts(r_1D,hist_r1D.binN);  % lin hist - autoscale bins to lim
hist_r1D.binCent=0.5*(hist_r1D.binEdge(1:end-1)+hist_r1D.binEdge(2:end));   % arith avg to bin center

% evaluate number density
nden_r1D=(hist_r1D.N)./(nShot*detQE*r_perp_area*diff(hist_r1D.binEdge));  % normalised for: shot, QE, phase space volume
nden_k1D=(hbar*tof/m_He)^3*nden_r1D;    % number density (real space) to far-field momentum density

if verbose>0    % plot
    % real-space density dist
    h_nr1D=figure();
    plot(1e3*hist_r1D.binCent,...
        1e-9*nden_r1D,'*-');    % scale units
    title('1D condensate number profile');
    xlabel('$r$ [mm]'); ylabel('$n(r,\overline{t})$ [mm$^{-3}$]');
    
    fname_str='nr1D_test';
    saveas(h_nr1D,[configs.files.dirout,fname_str,'.png']);
    saveas(h_nr1D,[configs.files.dirout,fname_str,'.fig']);
    
    % far-field momentum space
    h_nk1D=figure();
    plot(1e-6*r2k(hist_r1D.binCent),...
        1e18*nden_k1D,'*-');    % scale units
    title('1D condensate momentum profile');
    xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
    
    fname_str='nk1D';
    saveas(h_nk1D,[configs.files.dirout,fname_str,'.png']);
    saveas(h_nk1D,[configs.files.dirout,fname_str,'.fig']);
end


%%% LOG DIST
hist_lgk1D.binEdge=configs.hist.ed_lgk;     % get log-spaced edges
hist_lgk1D.binCent=sqrt(hist_lgk1D.binEdge(1:end-1).*hist_lgk1D.binEdge(2:end));   % GEOM avg bin centres
hist_lgk1D.N=histcounts(k_1D,hist_lgk1D.binEdge);	% use original data but bins are log-spaced

% evaluate number density
nden_lgk1D=(hist_lgk1D.N)./(nShot*detQE*k_perp_area*diff(hist_lgk1D.binEdge));  % normalised for: shot, QE, phase space volume

if verbose>0    % plot
    % far-field momentum space (log)
    h_nk1D_log=figure();
    loglog(1e-6*hist_lgk1D.binCent,...
        1e18*nden_lgk1D,'*-');     % scale units appropriately
    
    xlim([1e-1,2e1]);   %   limit x-axis to like Clement paper
    grid on;
    title('1D condensate momentum profile');
    xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
    
    fname_str='nk1D_loglog';
    saveas(h_nk1D_log,[configs.files.dirout,fname_str,'.png']);
    saveas(h_nk1D_log,[configs.files.dirout,fname_str,'.fig']);
end


%% n(k)k4 scaled plot
nk4_scaled=nden_lgk1D.*((hist_lgk1D.binCent).^4);

if verbose>0    % plot
    h_nk4=figure();
    semilogy(1e-6*hist_lgk1D.binCent,nk4_scaled,'*-');
    
    grid on;
    ylim([1e8,1e11]);       % y limits to like Clement PRL
    xlabel('$k$ [$\mu$m$^{-1}$]');
    ylabel('$k^{4}n_{\infty}(k)$ [m$^{-1}$]');
    
    fname_str='nk4_profile';
    saveas(h_nk4,[configs.files.dirout,fname_str,'.png']);
    saveas(h_nk4,[configs.files.dirout,fname_str,'.fig']);
end


%% Fit to density profile
% TODO


%% Save data
for i = 1:length(vars_save)
    if ~exist(vars_save{i},'var')
        warning(['Variable "',vars_save{i},'" does not exist.']);
        continue;
    end
    save(configs.files.saveddata,vars_save{i},'-append');     % to overcome version conflict
end