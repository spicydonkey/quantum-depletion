% Quantum depletion
% DKS 06/02/2017

clear all; close all; clc;

%%% USER INPUTS
path_config='C:\Users\HE BEC\Documents\MATLAB\quantum-depletion\config_test.m';

% vars to save to output
vars_save={'path_config',...
    'zxy_0','files_out',...
    'r_1D','k_1D',...
    'r_perp_area','k_perp_area'...
    'hist_r1D','nden_r1D','nden_k1D',...
    'hist_lgk1D','nden_lgk1D',...
    'nk4',...
    'nden_lgk_avg','nden_lgk_std','nden_lgk_se',...
    'k4_fit'...
    'hist_k_cyl_1D','nden_k_cyl_1D',...
    'nden_k_cyl_avg','nden_k_cyl_std','nden_k_cyl_se',...
    'k4cyl_fit'...
    };


%% Main
t_main_start=tic;

% load config
run(path_config);

% load misc params
hbar=configs.const.hbar;
m_He=configs.const.m_He;
tof=configs.const.tof;
vz=9.8*tof;     % atom free-fall vert v at detector hit for T-to-Z conversion;
detQE=configs.const.detect_qe;



hist_r1D.binN=configs.hist.nbin;    % get number of bins to set up auto

% Load TXY data and crop to region of interest
[txy_raw,files_out]=loadExpData(configs,verbose);
nShot_raw=size(txy_raw,1);


%% Pre-processing
% Build BEC locator box
bec_boxlim=cell(1,3);
for i=1:3
    bec_boxlim{i}=configs.bec.txy_pos(i)+configs.bec.box_fwidth(i)*[-0.5,0.5];
end

%%% Prepare complete set of oscillation cancelled ZXY data
% convert TXY to ZXY centred around condensate (oscillation compensation)
txy_0=cell(nShot_raw,1);        % oscillation compensated txy counts
zxy_0=cell(nShot_raw,1);        % T-Z conversion

% txy_cent=zeros(nShot_raw,3);    % centre of all captured TXY
bec_cent=zeros(nShot_raw,3);    % approx condensate centre from average of captured
num_txy_raw=zeros(1,nShot_raw); % number of counts in loaded shot

% build complete oscil cancelled data in the region of interest
for i=1:nShot_raw
%     txy_cent(i,:)=mean(txy_raw{i},1);       % approx condensate centre from average of captured
    % get mean position of counts captured in user specified box for
    % locating BEC
    [~,~,~,bec_cent(i,:)]=boxcull(txy_raw{i},bec_boxlim);   
    
    num_txy_raw(i)=size(txy_raw{i},1);      % number of counts in this shot
    txy_0{i}=txy_raw{i}-repmat(bec_cent(i,:),[num_txy_raw(i),1]);   % centre around self-average
    
    % T-Z conversion
    zxy_0{i}=txy_0{i};
    zxy_0{i}(:,1)=zxy_0{i}(:,1)*vz;     % TOF - dz at time of detection
end

%%% Plot far-field ZXY (summary)
nShotSumm=50;   % number of shots to plot as summary
if nShot_raw<nShotSumm
    nShotSumm=nShot_raw;
end
if verbose>1
    h_zxy_ff=figure();
    plot_zxy(zxy_0(1:nShotSumm),20,'k');    
    title('Condensate point cloud (Summary)');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    view(3);
    axis equal;
    
    % save plot
    fname_str='zxy_ff';
    saveas(h_zxy_ff,[configs.files.dirout,fname_str,'.png']);
    saveas(h_zxy_ff,[configs.files.dirout,fname_str,'.fig']);
end


%% METHOD 1: 1D-radial cylinder analysis
num_rot_angle=length(configs.axial_rot_angle);  % number of rotation angles to do analysis

Nmin_1D=configs.slice.mincount;     % minimum num reqd in 1D cyl to pass
N_slice_all=zeros(num_rot_angle,nShot_raw);     % store number captured in 1D slice

% phase space volume params
r_perp_area=pi*(configs.slice.cyl_rad^2);       % real area of integrated dims
k_perp_area=pi*(r2k(configs.slice.cyl_rad)^2);  % k area int area

% histogram params
hist_r1D.binN=configs.hist.nbin;    % get number of bins to set up (auto) for linear n profile

hist_lgk1D.binEdge=configs.hist.ed_lgk;     % get log-spaced edges
hist_lgk1D.binCent=sqrt(hist_lgk1D.binEdge(1:end-1).*hist_lgk1D.binEdge(2:end));   % GEOM avg bin centres
        
for i=1:num_rot_angle     % rotate whole zxy to sample 1D slice with angular shift
    x_thetha_tmp=configs.axial_rot_angle(i);    % rotation angle
    
    %% rotate zxy data around x
    zxy_xtheta_tmp=zxy_0;
    for j=1:nShot_raw
        z_tmp=zxy_xtheta_tmp{j}(:,1);
        y_tmp=zxy_xtheta_tmp{j}(:,3);
        
        zxy_xtheta_tmp{j}(:,1)=cos(x_thetha_tmp)*z_tmp-sin(x_thetha_tmp)*y_tmp;
        zxy_xtheta_tmp{j}(:,3)=sin(x_thetha_tmp)*z_tmp+cos(x_thetha_tmp)*y_tmp;
    end
    
    %% Take 1D slice in this configuration
    zxy_slice_tmp=cell(nShot_raw,1);    % 1D slice captured counts
    for j=1:nShot_raw
        [zxy_slice_tmp{j},~,n_captured]=cylindercull(zxy_xtheta_tmp{j},configs.slice.cyl_cent,...
            configs.slice.cyl_dim,configs.slice.cyl_orient);
        N_slice_all(i,j)=n_captured;         % save captured number
        
        % TODO - skipping low count in slice for now
        % could check for n_captured
        % filter shots based on counts in 1D slice
        % OLD CODE (VARS MAY BE OBSOLETE):
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if n_captured>Nmin_1D
%             txy_0{counter}=txy_0_temp;
%             zxy_0{counter}=zxy_0_temp;
%             zxy_slice{counter}=zxy_slice_temp;
%             fid_analysis(counter)=files_out.id_ok(i);
%             
%             counter=counter+1;
%         else
%             if verbose>0
%                 warning('1D filter: Low-count in file #%d. Discarding from further processing.',files_out.id_ok(i));
%             end
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         nShot=counter-1;    % number of shots passed from filtering for N in 1D slice
%         % clean up arrays
%         txy_0=txy_0(1:nShot);
%         zxy_0=zxy_0(1:nShot);
%         zxy_slice=zxy_slice(1:nShot);
%         fid_analysis=fid_analysis(1:nShot);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % Plot 1D slice point cloud
    if verbose>1
        if i==1
            h_zxy_slice=figure();
        end
        figure(h_zxy_slice); clf;
        plot_zxy(zxy_slice_tmp,100,'k');
        title(num2str(x_thetha_tmp));
        xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
        view(3);
        axis equal;
        
        % save plot
        fname_str=['zxy_slice_',num2str(i)];
        saveas(h_zxy_slice,[configs.files.dirout,fname_str,'.png']);
        saveas(h_zxy_slice,[configs.files.dirout,fname_str,'.fig']);
    end
    
    %% 1D cull and processing
    % collate and collapse all captured count to 1D
    r_1D=abs(vertcat(zxy_slice_tmp{:}));
    r_1D=r_1D(:,configs.slice.cyl_orient);	% 1D culled data in real-space

%     % collate all captured count in 1D-cyl slice from cond centre
%     r_1D=vertcat(zxy_slice_tmp{:});
%     % get magnitude to build array of 1D radial vectors
%     r_1D=sqrt(sum((r_1D.^2),2));    % NOTE: numerical error large for small r
    
    % r to k
    k_1D=r2k(r_1D);     % array of 1D k in [m^-1]
    
    %%%% Density profiling
    %% LINEAR DIST
    [hist_r1D.N{i},hist_r1D.binEdge{i}]=histcounts(r_1D,hist_r1D.binN);  % lin hist - autoscale bins to lim
    hist_r1D.binCent{i}=0.5*(hist_r1D.binEdge{i}(1:end-1)+hist_r1D.binEdge{i}(2:end));   % arith avg to bin center
    
    % evaluate number density
    nden_r1D{i}=(hist_r1D.N{i})./(nShot_raw*detQE*r_perp_area*diff(hist_r1D.binEdge{i}));  % normalised for: shot, QE, phase space volume
    nden_k1D{i}=(hbar*tof/m_He)^3*nden_r1D{i};    % number density (real space) to far-field momentum density
    
    if verbose>0    % plot
        % real-space density dist
        if i==1
            h_nr1D=figure();
            hold on; box on;
        end
        figure(h_nr1D);
        plot(1e3*hist_r1D.binCent{i},...
            1e-9*nden_r1D{i},'.-');    % scale units
        
        title('1D condensate number profile');
        xlabel('$r$ [mm]'); ylabel('$n(r,\overline{t})$ [mm$^{-3}$]');
        
        fname_str='nr1D_test';
        saveas(h_nr1D,[configs.files.dirout,fname_str,'.png']);
        saveas(h_nr1D,[configs.files.dirout,fname_str,'.fig']);
        
        % far-field momentum space
        if i==1
            h_nk1D=figure();
            hold on; box on;
        end
        figure(h_nk1D);
        plot(1e-6*r2k(hist_r1D.binCent{i}),...
            1e18*nden_k1D{i},'.-');    % scale units
        title('1D condensate momentum profile');
        xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
        
        if i==num_rot_angle     % save after all plotted
            fname_str='nk1D';
            saveas(h_nk1D,[configs.files.dirout,fname_str,'.png']);
            saveas(h_nk1D,[configs.files.dirout,fname_str,'.fig']);
        end
    end
    
    %% LOG DIST
    hist_lgk1D.N{i}=histcounts(k_1D,hist_lgk1D.binEdge);	% use original data but bins are log-spaced
    
    % evaluate number density
    nden_lgk1D{i}=(hist_lgk1D.N{i})./(nShot_raw*detQE*k_perp_area*diff(hist_lgk1D.binEdge));  % normalised for: shot, QE, phase space volume
    
    if verbose>0    % plot
        % far-field momentum space (log)
        if i==1
            h_nk1D_log=figure(); box on; hold on;
        end
        figure(h_nk1D_log);
        loglog(1e-6*hist_lgk1D.binCent,...
            1e18*nden_lgk1D{i},'.-');     % scale units appropriately
        
        xlim(1e6*configs.limit.k_com);
        ylim(1e18*configs.limit.kdensity);
        set(gca,'xScale','log');
        set(gca,'yScale','log');
        
        grid on;
        title('1D condensate momentum profile');
        xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
        
        if i==num_rot_angle     % save after all plotted
            fname_str='nk1D_loglog';
            saveas(h_nk1D_log,[configs.files.dirout,fname_str,'.png']);
            saveas(h_nk1D_log,[configs.files.dirout,fname_str,'.fig']);
        end
    end
    
    %% Evaluate n(k)k4 scaled plot
    nk4{i}=nden_lgk1D{i}.*((hist_lgk1D.binCent).^4);
    
    if verbose>0    % plot
        if i==1
            h_nk4=figure(); box on; hold on; 
        end
        figure(h_nk4);
        semilogy(1e-6*hist_lgk1D.binCent,nk4{i},'.-');
        
        ylim([1e8,1e11]);       % y limits to like Clement PRL
        set(gca,'yScale','log');
        
        grid on;
        xlabel('$k$ [$\mu$m$^{-1}$]');
        ylabel('$k^{4}n_{\infty}(k)$ [m$^{-1}$]');
        
        if i==num_rot_angle     % save after all plotted
            fname_str='nk4';
            saveas(h_nk4,[configs.files.dirout,fname_str,'.png']);
            saveas(h_nk4,[configs.files.dirout,fname_str,'.fig']);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% background k-density profile
bgd_r_perp_area=pi*(configs.bgd_cyl_rad^2);       % real area of integrated dims
bgd_k_perp_area=pi*(r2k(configs.bgd_cyl_rad)^2);  % k area int area

hist_r1D_bgd.binN=configs.hist.nbin;

hist_lgk1D_bgd.binEdge=configs.hist.ed_lgk;     % get log-spaced edges
hist_lgk1D_bgd.binCent=sqrt(hist_lgk1D_bgd.binEdge(1:end-1).*hist_lgk1D_bgd.binEdge(2:end));   % GEOM avg bin centres

bgd_zxy=cell(nShot_raw,1);
for i=1:nShot_raw
    bgd_zxy{i}=cylindercull(zxy_0{i},configs.bgd_cyl_cent,...
        configs.bgd_cyl_dim,configs.bgd_cyl_orient);
end
bgd_r_1D=abs(vertcat(bgd_zxy{:}));
bgd_r_1D=bgd_r_1D(:,configs.bgd_cyl_orient);	% 1D culled data in real-space
bgd_k_1D=r2k(bgd_r_1D);     % r to k

% Density profile
%%% Linear
[hist_r1D_bgd.N,hist_r1D_bgd.binEdge]=histcounts(bgd_r_1D,hist_r1D_bgd.binN);  % lin hist - autoscale bins to lim
hist_r1D_bgd.binCent=0.5*(hist_r1D_bgd.binEdge(1:end-1)+hist_r1D_bgd.binEdge(2:end));   % arith avg to bin center

nden_bgdr1D=(hist_r1D_bgd.N)./(nShot_raw*detQE*bgd_r_perp_area*diff(hist_r1D_bgd.binEdge));  % normalised for: shot, QE, phase space volume
nden_bgdk1D=(hbar*tof/m_He)^3*nden_bgdr1D;    % number density (real space) to far-field momentum density

% if verbose>0    % plot
%     % real-space density dist
%     figure(h_nr1D); hold on;
%     plot(1e3*hist_r1D_bgd.binCent,...
%         1e-9*nden_bgdr1D,'*-');    % scale units
%     title('1D condensate number profile');
%     xlabel('$r$ [mm]'); ylabel('$n(r,\overline{t})$ [mm$^{-3}$]');
%     
%     fname_str='nr1D_bgd';
%     saveas(h_nr1D,[configs.files.dirout,fname_str,'.png']);
%     saveas(h_nr1D,[configs.files.dirout,fname_str,'.fig']);
%     
%     % far-field momentum space
%     figure(h_nk1D); hold on;
%     plot(1e-6*r2k(hist_r1D_bgd.binCent),...
%         1e18*nden_bgdk1D,'*-');    % scale units
%     title('1D condensate momentum profile');
%     xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
%     
%     fname_str='nk1D_bgd';
%     saveas(h_nk1D,[configs.files.dirout,fname_str,'.png']);
%     saveas(h_nk1D,[configs.files.dirout,fname_str,'.fig']);
% end

%%% Log
hist_lgk1D_bgd.N=histcounts(bgd_k_1D,hist_lgk1D_bgd.binEdge);	% use original data but bins are log-spaced
nden_bgdlgk1D=(hist_lgk1D_bgd.N)./(nShot_raw*detQE*bgd_k_perp_area*diff(hist_lgk1D_bgd.binEdge));  % normalised for: shot, QE, phase space volume
n_bgd_smooth=smooth(nden_bgdlgk1D,ceil(length(hist_lgk1D_bgd.binCent)/10));

if verbose>0    % plot
    % far-field momentum space (log)
    figure(h_nk1D_log); hold on;
    loglog(1e-6*hist_lgk1D_bgd.binCent,...
        1e18*nden_bgdlgk1D,'.-');     % scale units appropriately
    
    grid on;
    title('1D condensate momentum profile');
    xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
    
    fname_str='nk1D_loglog_bgd';
    saveas(h_nk1D_log,[configs.files.dirout,fname_str,'.png']);
    saveas(h_nk1D_log,[configs.files.dirout,fname_str,'.fig']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% METHOD 2: Cylindrical sector for angular averaging
% get cylindrical sector params
cyl_dtheta=diff(configs.cylsect_theta_lims);
cyl_dktrans=2*r2k(configs.cylsect_trans_hwidth);

% convert BEC centered real-space counts to cylindrically symmetric coord system
% FORMAT: R_CYL = (rad_plane,theta,dist_transverse)
R0_cyl=cell(size(zxy_0));
for i=1:nShot_raw
    R0_cyl{i}=zeros(size(zxy_0{i}));
    R0_cyl{i}(:,1)=sqrt(sum(zxy_0{i}(:,[1,3]).^2,2));   % get in-plane radius [m]
    R0_cyl{i}(:,2)=atan2(zxy_0{i}(:,1),zxy_0{i}(:,3))+pi/2;     % theta origin is pointing "down" in Z [-pi,pi]
    R0_cyl{i}(:,3)=zxy_0{i}(:,2);       % transverse direction is in X [m]
end

% convert cylindrical R to k-space
k_cyl=R0_cyl;
for i=1:nShot_raw
    k_cyl{i}(:,1)=r2k(k_cyl{i}(:,1));
    k_cyl{i}(:,3)=r2k(k_cyl{i}(:,3));
end

% crop to region of interest for binning k
k_cyl_sect=cell(size(k_cyl));
for i=1:nShot_raw
    k_cyl_sect{i}=k_cyl{i};  % get everything
    
    % cull to angular lims
    theta_tmp=k_cyl_sect{i}(:,2);	% [-pi,pi]
    theta_tmp=wrapTo2Pi(theta_tmp-configs.cylsect_theta_lims(1));   % zero angle to start of lim
    ind_tmp=(theta_tmp<cyl_dtheta);   % indices lying in defined angular sector
    
    k_cyl_sect{i}=k_cyl_sect{i}(ind_tmp,:);   % cull
    
    % cull to transverse width
    k_trans_tmp=k_cyl_sect{i}(:,3);  % transverse k
    ind_tmp=(abs(k_trans_tmp)<r2k(configs.cylsect_trans_hwidth));   % indices lying in defined angular section
    
    k_cyl_sect{i}=k_cyl_sect{i}(ind_tmp,:);   % cull
end

% Get 1D-k
k_1D_cyl=cell(nShot_raw,1);    % get shot-wise cell array of k
for i=1:nShot_raw
    k_1D_cyl{i}=k_cyl_sect{i}(:,1);
end
% k_1D_cyl=vertcat(k_cyl_sect{:});
% k_1D_cyl=k_1D_cyl(:,1);

%%% Histogram
% get hist params
hist_k_cyl_1D.binEdge=configs.hist.ed_lgk;     % get log-spaced edges
hist_k_cyl_1D.binCent=sqrt(hist_k_cyl_1D.binEdge(1:end-1).*hist_k_cyl_1D.binEdge(2:end));   % GEOM avg bin centres

dk_cyl_volume=cyl_dktrans*cyl_dtheta*(hist_k_cyl_1D.binCent).*diff(hist_k_cyl_1D.binEdge);     % phase space volume in k

hist_k_cyl_1D.N=cell(nShot_raw,1);
nden_k_cyl_1D=cell(nShot_raw,1);
for i=1:nShot_raw
    hist_k_cyl_1D.N{i}=histcounts(k_1D_cyl{i},hist_k_cyl_1D.binEdge);   % single-shot
    nden_k_cyl_1D{i}=(hist_k_cyl_1D.N{i})./(detQE*dk_cyl_volume);  % normalised for: QE, phase space volume
end
% hist_k_cyl_1D.N=histcounts(k_1D_cyl,hist_k_cyl_1D.binEdge);

% evaluate number density
% nden_k_cyl_1D=(hist_k_cyl_1D.N)./(nShot_raw*detQE*dk_cyl_volume);  % normalised for: shot, QE, phase space volume

if verbose>0    % plot
    % far-field momentum space (log)
    h_nk_cyl_1D_log=figure();
    figure(h_nk_cyl_1D_log);
    hold on; box on;
    
    for i=1:nShot_raw
        plot(1e-6*hist_k_cyl_1D.binCent,...
            1e18*nden_k_cyl_1D{i},'.-');     % scale units appropriately
    end
    
    xlim(1e6*configs.limit.k_com);
    ylim(1e18*configs.limit.kdensity);
    set(gca,'xScale','log');
    set(gca,'yScale','log');
    
    grid on;
    title('1D condensate momentum profile - cylindrical capture');
    xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
    
    fname_str='nk_cyl_1D_loglog';
    saveas(h_nk_cyl_1D_log,[configs.files.dirout,fname_str,'.png']);
    saveas(h_nk_cyl_1D_log,[configs.files.dirout,fname_str,'.fig']);
end

% %%% Evaluate nk4
% nk4_cyl=nden_k_cyl_1D.*((hist_k_cyl_1D.binCent).^4);
% 
% if verbose>0    % plot
%     h_nk4_cyl=figure();
%     figure(h_nk4_cyl);
%     
%     semilogy(1e-6*hist_k_cyl_1D.binCent,nk4_cyl,'*-');
%     hold on;
%     
%     grid on;
%     ylim([1e8,1e11]);       % y limits to like Clement PRL
%     xlabel('$k$ [$\mu$m$^{-1}$]');
%     ylabel('$k^{4}n_{\infty}(k)$ [m$^{-1}$]');
%     
%     fname_str='nk4_cyl';
%     saveas(h_nk4_cyl,[configs.files.dirout,fname_str,'.png']);
%     saveas(h_nk4_cyl,[configs.files.dirout,fname_str,'.fig']);
% end


%% Statistical summary
% TODO - do for all other data
%%%% METHOD 1
%%% k distribution - log-log
nden_lgk_collated=vertcat(nden_lgk1D{:});   % collated k-density profile for all angles
nden_lgk_avg=mean(nden_lgk_collated,1);     % angular averaged farfield k profile
nden_lgk_std=std(nden_lgk_collated,1);      % standard deviation
nden_lgk_se=nden_lgk_std/sqrt(size(nden_lgk_collated,1));	% standard error

% Plot
h_nk_log_aa=figure(); 
hold on; box on;
mseb(1e-6*hist_lgk1D.binCent,1e18*nden_lgk_avg,...
    1e18*nden_lgk_std);  % NOTE: error in shaded error bar when error is larger than mean

% plot background count distribution
loglog(1e-6*hist_lgk1D_bgd.binCent,...
        1e18*n_bgd_smooth,'.-');     % scale units appropriately

xlim(1e6*configs.limit.k_com);
ylim(1e18*configs.limit.kdensity);
set(gca,'xScale','log');    % loglog scale
set(gca,'yScale','log');

grid on;
title('Angular averaged - 1D k profile');
xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');

fname_str='nk1D_loglog_aa';
saveas(h_nk_log_aa,[configs.files.dirout,fname_str,'.png']);
saveas(h_nk_log_aa,[configs.files.dirout,fname_str,'.fig']);


%%% n(k)k4
nk4_collated=vertcat(nk4{:});       % collated nk4 profile for all angles
nk4_avg=mean(nk4_collated,1);       % angular averaged
nk4_std=std(nk4_collated,1);
nk4_se=nk4_std/sqrt(size(nk4_collated,1));  % standard error

% Plot
h_nk4_aa=figure();
hold on; box on;
mseb(1e-6*hist_lgk1D.binCent,nk4_avg,...
    nk4_std);        % NOTE: error in shaded error bar when error is larger than mean

set(gca,'yScale','log');    % y scale log

grid on;
ylim([1e8,1e11]);       % y limits to like Clement PRL
xlabel('$k$ [$\mu$m$^{-1}$]');
ylabel('$k^{4}n_{\infty}(k)$ [m$^{-1}$]');

fname_str='nk4_aa';
saveas(h_nk4_aa,[configs.files.dirout,fname_str,'.png']);
saveas(h_nk4_aa,[configs.files.dirout,fname_str,'.fig']);

%%%% METHOD 2
nden_k_cyl_collated=vertcat(nden_k_cyl_1D{:});   % collated k-density profile for all angles
nden_k_cyl_avg=mean(nden_k_cyl_collated,1);     % angular averaged farfield k profile
nden_k_cyl_std=std(nden_k_cyl_collated,1);      % standard deviation
nden_k_cyl_se=nden_k_cyl_std/sqrt(size(nden_k_cyl_collated,1));	% standard error

% Plot
h_nk_cyl_avg=figure();
hold on; box on;
mseb(1e-6*hist_k_cyl_1D.binCent,1e18*nden_k_cyl_avg,...
    1e18*nden_k_cyl_std);  % NOTE: error in shaded error bar when error is larger than mean

% plot smoothed background count distribution
loglog(1e-6*hist_lgk1D_bgd.binCent,...
        1e18*n_bgd_smooth,'.-');     % scale units appropriately
    
xlim(1e6*configs.limit.k_com);
ylim(1e18*configs.limit.kdensity);
set(gca,'xScale','log');    % loglog scale
set(gca,'yScale','log');

grid on;
title('Angular averaged - 1D k profile');
xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');

fname_str='nk1D_cyl_avg';
saveas(h_nk_cyl_avg,[configs.files.dirout,fname_str,'.png']);
saveas(h_nk_cyl_avg,[configs.files.dirout,fname_str,'.fig']);


%% Fit density profile
%%% METHOD 1
% Negative power law fit to k-distribution at large k
[~,I_qd]=min(abs(hist_lgk1D.binCent-configs.fit.k_min));  % get index from which to fit QD neg-power law
k4_fit.QD.k=hist_lgk1D.binCent(I_qd:end);   % store data used for fitting
k4_fit.QD.nk=nden_lgk_avg(I_qd:end);

% Call the Matlab fitting routine
k4_fit.QD.fit=fitnlm(k4_fit.QD.k,k4_fit.QD.nk,...
    configs.fit.fun_negpowk,configs.fit.param0,...
    'CoefficientNames',configs.fit.fun_coefname,...
    'Options',configs.fit.opt);

% Summarise fit
disp(k4_fit.QD.fit);

% build a sample of the fitted model's profile
ratio_extrap=1.5;
k4_fit.QD.k_fit=logspace(log10(min(k4_fit.QD.k)/ratio_extrap),log10(ratio_extrap*max(k4_fit.QD.k)),1000);   % indep var to evaluate fitted function
k4_fit.QD.nk_fit=feval(k4_fit.QD.fit,k4_fit.QD.k_fit);  % evaluate fitted model

% Plot
figure(h_nk_log_aa); hold on;
loglog(1e-6*k4_fit.QD.k_fit,1e18*k4_fit.QD.nk_fit,'k--');

% Save plot
fname_str='nk_fit';
saveas(h_nk_log_aa,[configs.files.dirout,fname_str,'.png']);
saveas(h_nk_log_aa,[configs.files.dirout,fname_str,'.fig']);


%%% METHOD 2
% I_qd is shared
% [~,I_qd]=min(abs(hist_k_cyl_1D.binCent-configs.fit.k_min));  % get index from which to fit QD neg-power law
k4cyl_fit.QD.k=hist_k_cyl_1D.binCent(I_qd:end);   % store data used for fitting
k4cyl_fit.QD.nk=nden_k_cyl_avg(I_qd:end);

k4cyl_fit.QD.fit=fitnlm(k4cyl_fit.QD.k,k4cyl_fit.QD.nk,...
    configs.fit.fun_negpowk,configs.fit.param0,...
    'CoefficientNames',configs.fit.fun_coefname,...
    'Options',configs.fit.opt);

% Summarise fit
disp(k4cyl_fit.QD.fit);

% build a sample of the fitted model's profile
% ratio_extrap=1.5;
k4cyl_fit.QD.k_fit=logspace(log10(min(k4cyl_fit.QD.k)/ratio_extrap),log10(ratio_extrap*max(k4cyl_fit.QD.k)),1000);   % indep var to evaluate fitted function
k4cyl_fit.QD.nk_fit=feval(k4cyl_fit.QD.fit,k4cyl_fit.QD.k_fit);  % evaluate fitted model

% Plot
figure(h_nk_cyl_avg); hold on;
loglog(1e-6*k4cyl_fit.QD.k_fit,1e18*k4cyl_fit.QD.nk_fit,'k--');

% Save plot
fname_str='nk_cyl_fit';
saveas(h_nk_cyl_avg,[configs.files.dirout,fname_str,'.png']);
saveas(h_nk_cyl_avg,[configs.files.dirout,fname_str,'.fig']);


%% Save data
for i = 1:length(vars_save)
    if ~exist(vars_save{i},'var')
        warning(['Variable "',vars_save{i},'" does not exist.']);
        continue;
    end
    save(configs.files.saveddata,vars_save{i},'-append');     % to overcome version conflict
end

%% END
t_main_end=toc(t_main_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
disp('===================ALL TASKS COMPLETED===================');