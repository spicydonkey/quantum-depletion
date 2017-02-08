% Quantum depletion
% DKS 06/02/2017

clear all; close all; clc;

%%% USER INPUTS
path_config='C:\Users\HE BEC\Documents\MATLAB\quantum-depletion\config_060217_test.m';

% vars to save to output
vars_save={'path_config',...
    'zxy_0','files_out',...
    'r_1D','k_1D',...
    'r_perp_area','k_perp_area'...
    'hist_r1D','nden_r1D','nden_k1D',...
    'hist_lgk1D','nden_lgk1D',...
    'nk4',...
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

%%% evaluate params for analysis
r_perp_area=pi*(configs.slice.cyl_rad^2);       % real area of integrated dims
k_perp_area=pi*(r2k(configs.slice.cyl_rad)^2);  % k area int area

hist_r1D.binN=configs.hist.nbin;    % get number of bins to set up auto

% Load TXY data and crop to region of interest
[txy_raw,files_out]=loadExpData(configs,verbose);
nShot_raw=size(txy_raw,1);


%% Pre-processing
%%% Prepare complete set of oscillation cancelled ZXY data
% convert TXY to ZXY centred around condensate (oscillation compensation)
txy_0=cell(nShot_raw,1);        % oscillation compensated txy counts
zxy_0=cell(nShot_raw,1);        % T-Z conversion

txy_cent=zeros(nShot_raw,3);    % approx condensate centre from average of captured
num_txy_raw=zeros(1,nShot_raw); % number of counts in loaded shot

% build complete oscil cancelled data in the region of interest
for i=1:nShot_raw
    txy_cent(i,:)=mean(txy_raw{i},1);       % approx condensate centre from average of captured
    num_txy_raw(i)=size(txy_raw{i},1);      % number of counts in this shot
    txy_0{i}=txy_raw{i}-repmat(txy_cent(i,:),[num_txy_raw(i),1]);   % centre around self-average
    
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


%% 1D calculation
num_rot_angle=length(configs.axial_rot_angle);

Nmin_1D=configs.slice.mincount;     % minimum count in 1D slice to pass
N_slice_all=zeros(num_rot_angle,nShot_raw);        % number captured in 1D slice

hist_r1D.binN=configs.hist.nbin;    % get number of bins to set up auto

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
        % PREV CODE:
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
            hold on;
        end
        figure(h_nr1D);
        plot(1e3*hist_r1D.binCent{i},...
            1e-9*nden_r1D{i},'*-');    % scale units
        title('1D condensate number profile');
        xlabel('$r$ [mm]'); ylabel('$n(r,\overline{t})$ [mm$^{-3}$]');
        
        fname_str='nr1D_test';
        saveas(h_nr1D,[configs.files.dirout,fname_str,'.png']);
        saveas(h_nr1D,[configs.files.dirout,fname_str,'.fig']);
        
        % far-field momentum space
        if i==1
            h_nk1D=figure();
            hold on;
        end
        figure(h_nk1D);
        plot(1e-6*r2k(hist_r1D.binCent{i}),...
            1e18*nden_k1D{i},'*-');    % scale units
        title('1D condensate momentum profile');
        xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
        
        fname_str='nk1D';
        saveas(h_nk1D,[configs.files.dirout,fname_str,'.png']);
        saveas(h_nk1D,[configs.files.dirout,fname_str,'.fig']);
    end
    
    %% LOG DIST
    hist_lgk1D.N{i}=histcounts(k_1D,hist_lgk1D.binEdge);	% use original data but bins are log-spaced
    
    % evaluate number density
    nden_lgk1D{i}=(hist_lgk1D.N{i})./(nShot_raw*detQE*k_perp_area*diff(hist_lgk1D.binEdge));  % normalised for: shot, QE, phase space volume
    
    if verbose>0    % plot
        % far-field momentum space (log)
        if i==1
            h_nk1D_log=figure();
        end
        figure(h_nk1D_log);
        loglog(1e-6*hist_lgk1D.binCent,...
            1e18*nden_lgk1D{i},'*-');     % scale units appropriately
        hold on;
        
        xlim([1e-1,2e1]);   %   limit x-axis to like Clement paper
        grid on;
        title('1D condensate momentum profile');
        xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
        
        fname_str='nk1D_loglog';
        saveas(h_nk1D_log,[configs.files.dirout,fname_str,'.png']);
        saveas(h_nk1D_log,[configs.files.dirout,fname_str,'.fig']);
    end
    
    %% Evaluate n(k)k4 scaled plot
    nk4{i}=nden_lgk1D{i}.*((hist_lgk1D.binCent).^4);
    
    if verbose>0    % plot
        if i==1
            h_nk4=figure();
        end
        figure(h_nk4);
        semilogy(1e-6*hist_lgk1D.binCent,nk4{i},'*-');
        hold on;
        
        grid on;
        ylim([1e8,1e11]);       % y limits to like Clement PRL
        xlabel('$k$ [$\mu$m$^{-1}$]');
        ylabel('$k^{4}n_{\infty}(k)$ [m$^{-1}$]');
        
        fname_str='nk4';
        saveas(h_nk4,[configs.files.dirout,fname_str,'.png']);
        saveas(h_nk4,[configs.files.dirout,fname_str,'.fig']);
    end
end


%% Statistical summary
% TODO - do for all other data
%%% k distribution - log-log
nden_lgk_collated=vertcat(nden_lgk1D{:});   % collated k-density profile for all angles
nden_lgk_avg=mean(nden_lgk_collated,1);     % angular averaged farfield k profile
nden_lgk_se=std(nden_lgk_collated,1)/sqrt(size(nden_lgk_collated,1));   % standard error

% Plot
h_nk_log_aa=figure();
loglog(1e-6*hist_lgk1D.binCent,...
    1e18*nden_lgk_avg,'o');     % scale units appropriately

xlim([1e-1,2e1]);   %   limit x-axis to like Clement paper
grid on;
title('Angular averaged - 1D k profile');
xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');

fname_str='nk1D_loglog_aa';
saveas(h_nk_log_aa,[configs.files.dirout,fname_str,'.png']);
saveas(h_nk_log_aa,[configs.files.dirout,fname_str,'.fig']);


%%% n(k)k4
nk4_collated=vertcat(nk4{:});       % collated nk4 profile for all angles
nk4_avg=mean(nk4_collated,1);       % angular averaged
nk4_se=std(nk4_collated,1)/sqrt(size(nk4_collated,1));  % standard error

% Plot
h_nk4_aa=figure();
semilogy(1e-6*hist_lgk1D.binCent,nk4_avg,'o');

grid on;
ylim([1e8,1e11]);       % y limits to like Clement PRL
xlabel('$k$ [$\mu$m$^{-1}$]');
ylabel('$k^{4}n_{\infty}(k)$ [m$^{-1}$]');

fname_str='nk4_aa';
saveas(h_nk4_aa,[configs.files.dirout,fname_str,'.png']);
saveas(h_nk4_aa,[configs.files.dirout,fname_str,'.fig']);

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