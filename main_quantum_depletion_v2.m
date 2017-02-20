% Quantum depletion
% DKS 06/02/2017

clear all; close all; clc;

%%% USER INPUTS
path_config='C:\Users\HE BEC\Documents\MATLAB\quantum-depletion\config_v2.m';
path_param_log='C:\Users\HE BEC\Documents\lab\quantum-depletion\exp4\log_test.txt';

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
nShot_all=size(txy_raw,1);

% Load log file with param settings ID
% table is row vectors of [file id, param set id]
% TODO: make modular to accept datetime format, comment, etc
% table_param_id=csvread(path_param_log);
% TODO: fix this HACK - assumes count number comparison to distinguish
table_param_id=zeros(nShot_all,2);
table_param_id(:,1)=files_out.id_ok;
for idx=1:nShot_all
    if size(txy_raw{idx},1)>1000
        table_param_id(idx,2)=1;
    else
        table_param_id(idx,2)=2;
    end
end

%% Pre-processing
txy_cat=cell(configs.paramset,1);     % TXY data categorised
settings_idx=cell(configs.paramset,1);

% categorise TXY files
for idx=1:configs.paramset
    settings_idx{idx}=find(table_param_id(:,2)==idx);
%     txy_0{idx}=cell(sum(table_param_id(:,2)==idx),1);
%     txy_0{idx}=cell(numel(settings_idx{idx}),1);
    txy_cat{idx}=txy_raw(settings_idx{idx});  % populate this param set
end

% %% Categorize Raw TXY data
% txy_raw_cat=cell(configs.paramset,1);
% % TODO: make general
% counter1=1;
% counter2=1;
% for i=1:size(txy_raw,1)
%     if mod(files_out.id_ok(i),2)==1
%         % Cat 1
%         txy_raw_cat{1}{counter1,1}=txy_raw{i};
%         counter1=counter1+1;
%     else
%         % Cat 2
%         txy_raw_cat{2}{counter2,1}=txy_raw{i};
%         counter2=counter2+1;
%     end
% end

% %% TODO: QUICK HACK FOR NOW
% % update all new data into old var settings
% txy_raw=txy_raw_cat{1};     % do same analysis but use only the param set #1
% 
% % re-evaluate
% nShot_all=size(txy_raw,1);


%%% Prepare complete set of oscillation cancelled ZXY data
% convert TXY to ZXY centred around condensate (oscillation compensation)
txy_0=cell(configs.paramset,1);        % oscillation compensated txy counts
zxy_0=cell(configs.paramset,1);        % T-Z conversion

for paridx=1:configs.paramset
    txy_0{paridx}=cell(size(txy_cat{paridx}));    % initialise cell array size
    zxy_0{paridx}=txy_0{paridx};
end

%% Param set #1
txy_temp=txy_cat{1};      % the signal BEC set
nShot_temp=numel(txy_temp);     % number of shots in set

% Build BEC locator box
bec_boxlim=cell(1,3);
for i=1:3
    bec_boxlim{i}=configs.bec.txy_pos(i)+configs.bec.box_fwidth(i)*[-0.5,0.5];
end

% txy_cent=zeros(nShot_raw,3);    % centre of all captured TXY
bec_cent=zeros(nShot_temp,3);    % approx condensate centre from average of captured
num_txy_raw=zeros(1,nShot_temp); % number of counts in loaded shot

% build complete oscil cancelled data in the region of interest
for i=1:nShot_temp
%     txy_cent(i,:)=mean(txy_raw{i},1);       % approx condensate centre from average of captured
    % get mean position of counts captured in user specified box for
    % locating BEC
    [~,~,~,bec_cent(i,:)]=boxcull(txy_temp{i},bec_boxlim);
    
    num_txy_raw(i)=size(txy_temp{i},1);      % number of counts in this shot
    txy_0{1}{i}=txy_temp{i}-repmat(bec_cent(i,:),[num_txy_raw(i),1]);   % centre around self-average
    
    % T-Z conversion
    zxy_0{1}{i}=txy_0{1}{i};
    zxy_0{1}{i}(:,1)=zxy_0{1}{i}(:,1)*vz;     % TOF - dz at time of detection
end

% get mean absolute position of BEC
% TODO - maybe cull bad shots
avg_bec_cent=nanmean(bec_cent,1);
std_bec_cent=nanstd(bec_cent,1);


%% Param set #2
txy_temp=txy_cat{2};    % the background counts set (RF sweep off)
nShot_temp=numel(txy_temp);     % number of shots in set
num_txy_raw=zeros(1,nShot_temp); % number of counts in loaded shot

% Centre to mean BEC position
for i=1:nShot_temp
    num_txy_raw(i)=size(txy_temp{i},1);      % number of counts in this shot
    txy_0{2}{i}=txy_temp{i}-repmat(avg_bec_cent,[num_txy_raw(i),1]);   % centre around BEC avg
    
    % T-Z conversion
    zxy_0{2}{i}=txy_0{2}{i};
    zxy_0{2}{i}(:,1)=zxy_0{2}{i}(:,1)*vz;     % TOF - dz at time of detection
end


%% Plot far-field ZXY (summary)
nShotSumm=50;   % number of shots to plot as summary
if min(size(zxy_0{1},1),size(zxy_0{2},1))<nShotSumm
    nShotSumm=min(size(zxy_0{1},1),size(zxy_0{2},1));
end
if verbose>1
    h_zxy_ff=figure();
    
    subplot(1,2,1);
    plot_zxy(zxy_0{1}(1:nShotSumm),20,'k');    
    title('Condensate point cloud (Summary)');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    view(3);
    axis equal;
    
    subplot(1,2,2);
    plot_zxy(zxy_0{2}(1:nShotSumm),20,'k');    
    title('Condensate point cloud (Summary)');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    view(3);
    axis equal;
    
    % save plot
    fname_str='zxy_ff';
    saveas(h_zxy_ff,[configs.files.dirout,fname_str,'.png']);
    saveas(h_zxy_ff,[configs.files.dirout,fname_str,'.fig']);
end


%% Cylindrical sector for angular averaging
% get cylindrical sector params
cyl_dtheta=diff(configs.cylsect_theta_lims);
cyl_dktrans=2*r2k(configs.cylsect_trans_hwidth);

for idxparam=1:configs.paramset
    
    % convert BEC centered real-space counts to cylindrically symmetric coord system
    % FORMAT: R_CYL = (rad_plane,theta,dist_transverse)
    R0_cyl{idxparam}=cell(size(zxy_0{idxparam}));
    nShot_this=size(zxy_0{idxparam},1);
    for i=1:nShot_this
        R0_cyl{idxparam}{i}=zeros(size(zxy_0{idxparam}{i}));
        R0_cyl{idxparam}{i}(:,1)=sqrt(sum(zxy_0{idxparam}{i}(:,[1,3]).^2,2));   % get in-plane radius [m]
        R0_cyl{idxparam}{i}(:,2)=atan2(zxy_0{idxparam}{i}(:,1),zxy_0{idxparam}{i}(:,3))+pi/2;     % theta origin is pointing "down" in Z [-pi,pi]
        R0_cyl{idxparam}{i}(:,3)=zxy_0{idxparam}{i}(:,2);       % transverse direction is in X [m]
    end
    
    % convert cylindrical R to k-space
    k_cyl{idxparam}=R0_cyl{idxparam};
    for i=1:nShot_this
        k_cyl{idxparam}{i}(:,1)=r2k(k_cyl{idxparam}{i}(:,1));
        k_cyl{idxparam}{i}(:,3)=r2k(k_cyl{idxparam}{i}(:,3));
    end
    
    % crop to region of interest for binning k
    k_cyl_sect{idxparam}=cell(size(k_cyl{idxparam}));
    for i=1:nShot_this
        k_cyl_sect{idxparam}{i}=k_cyl{idxparam}{i};  % get everything
        
        % cull to angular lims
        theta_tmp=k_cyl_sect{idxparam}{i}(:,2);	% [-pi,pi]
        theta_tmp=wrapTo2Pi(theta_tmp-configs.cylsect_theta_lims(1));   % zero angle to start of lim
        ind_tmp=(theta_tmp<cyl_dtheta);   % indices lying in defined angular sector
        
        k_cyl_sect{idxparam}{i}=k_cyl_sect{idxparam}{i}(ind_tmp,:);   % cull
        
        % cull to transverse width
        k_trans_tmp=k_cyl_sect{idxparam}{i}(:,3);  % transverse k
        ind_tmp=(abs(k_trans_tmp)<r2k(configs.cylsect_trans_hwidth));   % indices lying in defined angular section
        
        k_cyl_sect{idxparam}{i}=k_cyl_sect{idxparam}{i}(ind_tmp,:);   % cull
    end
    
    % Get 1D-k
    k_1D_cyl{idxparam}=cell(nShot_this,1);    % get shot-wise cell array of k
    for i=1:nShot_this
        k_1D_cyl{idxparam}{i}=k_cyl_sect{idxparam}{i}(:,1);
    end
    % k_1D_cyl=vertcat(k_cyl_sect{:});
    % k_1D_cyl=k_1D_cyl(:,1);
    
    %COLLATE!
    k_1D_cyl{idxparam}=vertcat(k_1D_cyl{idxparam}{:});
    
    %%% Histogram
    % get hist params
    hist_k_cyl_1D.binEdge{idxparam}=configs.hist.ed_lgk;     % get log-spaced edges
    hist_k_cyl_1D.binCent{idxparam}=sqrt(hist_k_cyl_1D.binEdge{idxparam}(1:end-1).*hist_k_cyl_1D.binEdge{idxparam}(2:end));   % GEOM avg bin centres
    
    dk_cyl_volume{idxparam}=cyl_dktrans*cyl_dtheta*(hist_k_cyl_1D.binCent{idxparam}).*diff(hist_k_cyl_1D.binEdge{idxparam});     % phase space volume in k
    
    % hist_k_cyl_1D.N=cell(nShot_raw,1);
    % nden_k_cyl_1D=cell(nShot_raw,1);
    % for i=1:nShot_raw
    %     hist_k_cyl_1D.N{i}=histcounts(k_1D_cyl{i},hist_k_cyl_1D.binEdge);   % single-shot
    %     nden_k_cyl_1D{i}=(hist_k_cyl_1D.N{i})./(detQE*dk_cyl_volume);  % normalised for: QE, phase space volume
    % end
    hist_k_cyl_1D.N{idxparam}=histcounts(k_1D_cyl{idxparam},hist_k_cyl_1D.binEdge{idxparam});     % collated
    
    % evaluate number density
    nden_k_cyl_1D{idxparam}=(hist_k_cyl_1D.N{idxparam})./(nShot_this*detQE*dk_cyl_volume{idxparam});  % normalised for: shot, QE, phase space volume
end

% subtract background
nk_bgd_free=nden_k_cyl_1D{1}-nden_k_cyl_1D{2};
nk_bgd_free(nk_bgd_free<=0)=NaN;    % handle negative density

% clean data
idx_isnum=~isnan(nk_bgd_free);
nk_bgd_free_clean=nk_bgd_free(idx_isnum);
k_clean=hist_k_cyl_1D.binCent{1}(idx_isnum);    % binCent needs to be common

if verbose>0    % plot
    % far-field momentum space (log)
    h_nk_cyl_1D_log=figure();
    figure(h_nk_cyl_1D_log);
    hold on; box on;
    
    for idxparam=1:configs.paramset
        plot(1e-6*hist_k_cyl_1D.binCent{idxparam},...
            1e18*nden_k_cyl_1D{idxparam},'.-');     % scale units appropriately
        
        xlim(1e6*configs.limit.k_com);
        ylim(1e18*configs.limit.kdensity);
        set(gca,'xScale','log');
        set(gca,'yScale','log');
        
        grid on;
        title('1D condensate momentum profile - cylindrical capture');
        xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
    end
    
    % plot background subtracted density profile
    plot(1e-6*hist_k_cyl_1D.binCent{idxparam},...
            1e18*nk_bgd_free,'-');     % scale units appropriately
    
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


%% Smooth profile
nk_sm_avg=smooth(nk_bgd_free_clean,configs.smooth.nspan,'moving');  % simple moving average filter

% Gaussian smoothing
g_filt=gausswin(configs.smooth.nspan);
g_filt=g_filt/sum(g_filt);
nk_sm_gauss=conv(nk_bgd_free_clean,g_filt,'same');

nk_sm_std=movingstd(nk_bgd_free_clean,configs.smooth.nspan,'central');   % moving std

% smoothing in logspace


% plot smoothed data
if verbose>0    % plot
    h_nk_sm=figure();
    figure(h_nk_sm);
    hold on; box on;
    
    shadedErrorBar(1e-6*k_clean,1e18*nk_sm_avg,1e18*nk_sm_std,'k');
    
    shadedErrorBar(1e-6*k_clean,1e18*nk_sm_gauss,1e18*nk_sm_std,'b');
    
    set(gca,'xScale','log');
    set(gca,'yScale','log');
    axis auto;
    grid on;
    
    title('smoothed momentum density profile');
    xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
    
    fname_str='nk_smoothed';
    saveas(h_nk_sm,[configs.files.dirout,fname_str,'.png']);
    saveas(h_nk_sm,[configs.files.dirout,fname_str,'.fig']);
end


%% Fit density profile
% common
ratio_extrap=1;

% fit to raw data
[~,I_qd]=min(abs(hist_k_cyl_1D.binCent{idxparam}-configs.fit.k_min));  % get index from which to fit QD neg-power law

% get fitting data
k4cyl_fit.QD.k_log=log(hist_k_cyl_1D.binCent{idxparam}(I_qd:end));   
k4cyl_fit.QD.nk_log=real(log(nk_bgd_free(I_qd:end)));
k4cyl_fit.QD.nk_log(isinf(k4cyl_fit.QD.nk_log))=NaN;    % infinity to NaN

k4cyl_fit.QD.fit=fitnlm(k4cyl_fit.QD.k_log,k4cyl_fit.QD.nk_log,...
    configs.fit.fun_negpowk,configs.fit.param0,...
    'CoefficientNames',configs.fit.fun_coefname,...
    'Options',configs.fit.opt);

% Summarise fit
disp(k4cyl_fit.QD.fit);

% build a sample of the fitted model's profile
% ratio_extrap=1.5;
k4cyl_fit.QD.k_log_fit=linspace((min(k4cyl_fit.QD.k_log)),max(k4cyl_fit.QD.k_log),1000);   % indep var to evaluate fitted function
k4cyl_fit.QD.nk_log_fit=feval(k4cyl_fit.QD.fit,k4cyl_fit.QD.k_log_fit);  % evaluate fitted model

% Plot
figure(h_nk_cyl_1D_log); hold on;
plot(1e-6*exp(k4cyl_fit.QD.k_log_fit),1e18*exp(k4cyl_fit.QD.nk_log_fit),'k--');

% auto axis
axis auto;

% Save plot
fname_str='nk_cyl_fit';
% saveas(h_nk_cyl_avg,[configs.files.dirout,fname_str,'.png']);
% saveas(h_nk_cyl_avg,[configs.files.dirout,fname_str,'.fig']);

saveas(h_nk_cyl_1D_log,[configs.files.dirout,fname_str,'.png']);
saveas(h_nk_cyl_1D_log,[configs.files.dirout,fname_str,'.fig']);


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