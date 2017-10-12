% Quantum depletion
% DKS 06/02/2017

clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUTS
fname_config='config_run6.m';       % name of config file to use (located in $QD_REPO/configs)

% path_config='C:\Users\HE BEC\Documents\MATLAB\quantum-depletion\configs\config_run7.m';

% NOTE: getting param id from logfile is not implemented yet
% path_param_log='C:\Users\HE BEC\Documents\lab\quantum-depletion\exp5\log_test.txt';   

%%% END OF USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% build paths
thispath=fileparts(mfilename('fullpath'));      % get directory to current codebase
path_config=fullfile(thispath,'configs',fname_config);  % build path to config file

% vars to save to output
vars_save={'configs','path_config',...
    'zxy_0','files_out',...
    'k_1D_cyl',...
    'hist_k_cyl_1D','nden_k_cyl_1D',...
    'nk','nk_se'...
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
    txy_cat{idx}=txy_raw(settings_idx{idx});  % populate this param set
end

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
    % get mean position of counts captured in user specified box for
    % locating BEC
    % TODO: strongly saturated BEC will be blobby and averaging will be
    % biased
    [~,~,~,bec_cent(i,:)]=boxcull(txy_temp{i},bec_boxlim);
    
    num_txy_raw(i)=size(txy_temp{i},1);      % number of counts in this shot
    txy_0{1}{i}=txy_temp{i}-repmat(bec_cent(i,:),[num_txy_raw(i),1]);   % centre around self-average
    
    % T-Z conversion
    zxy_0{1}{i}=txy_0{1}{i};
    zxy_0{1}{i}(:,1)=zxy_0{1}{i}(:,1)*vz;     % TOF - dz at time of detection
end

% get mean absolute position of BEC
% TODO - how can there may be bad shots that produce NaN as bec_cent?
avg_bec_cent=mean(bec_cent,1);
std_bec_cent=std(bec_cent,1);


%% Param set #2
if configs.paramset>1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEBUG FOR FLAT BACKGROUND
    if do_flat
        run('flat_background.m');
    end
    if exist('TXY_FLAT','var')
        txy_temp=TXY_FLAT;
        avg_bec_cent=[0,0,0];
    else
        txy_temp=txy_cat{2};    % the background counts set (RF sweep off)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
end


%% Plot far-field ZXY (summary)
% TODO: handle multiple param sets
if verbose>1
%     nShotSumm=50;   % number of shots to plot as summary
%     if min(size(zxy_0{1},1),size(zxy_0{2},1))<nShotSumm
%         nShotSumm=min(size(zxy_0{1},1),size(zxy_0{2},1));
%     end
    
    h_zxy_ff=figure();
    figure(h_zxy_ff);
%     subplot(1,2,1);
%     plot_zxy(zxy_0{1}(1:nShotSumm),20,'k');    
%     title('Condensate point cloud (Summary)');
%     xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
%     view(3);
%     axis equal;
%     
%     subplot(1,2,2);
%     plot_zxy(zxy_0{2}(1:nShotSumm),20,'k');    
%     title('Condensate point cloud (Summary)');
%     xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
%     view(3);
%     axis equal;
    
    plot_zxy(zxy_0{1},configs.num_count_disp,20,'k');
    
%     title_str=sprintf('');
%     title(title_str);
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    view(3);
    axis equal;
    
%     % save plot
%     fname_str='zxy_ff';
%     saveas(h_zxy_ff,[configs.files.dirout,fname_str,'.png']);
%     saveas(h_zxy_ff,[configs.files.dirout,fname_str,'.fig']);
end

%% Far-field k-vectors
% TODO: can use k_zxy as data to update following procedures

% Build data in cartesian k-vector array
k_zxy=cell(configs.paramset,1);     % data in Cartesian k-vector array
for idxparam=1:configs.paramset
    nShot_this=size(zxy_0{idxparam},1);
    k_zxy{idxparam}=cell(nShot_this,1);
    for iShot=1:nShot_this
        k_zxy{idxparam}{iShot}=r2k(zxy_0{idxparam}{iShot});
    end
end

%% Plot Far-field k-vectors (summary)
% TODO: handle multiple param sets
if verbose>1    
    h_k_ff=figure();
    figure(h_k_ff);
    
    plot_zxy(k_zxy{1},configs.num_count_disp,20,'k');
    
    xlabel('$k_X$ [m$^{-1}$]'); ylabel('$k_Y$ [m$^{-1}$]'); zlabel('$k_Z$ [m$^{-1}$]');
    view(3);
    axis equal;
    
%     % save plot
%     fname_str='k_ff';
%     saveas(h_k_ff,[configs.files.dirout,fname_str,'.png']);
%     saveas(h_k_ff,[configs.files.dirout,fname_str,'.fig']);
end


%% Cart-Cyl(axis=X) coord transform
% convert BEC centered real-space counts to cylindrically symmetric coord system
% FORMAT: R_CYL = (rad_plane,theta,dist_transverse)
R0_cyl=cell(configs.paramset,1);    % init cell array
k_cyl=cell(configs.paramset,1);

for idxparam=1:configs.paramset
    R0_cyl{idxparam}=cell(size(zxy_0{idxparam}));
    nShot_this=size(zxy_0{idxparam},1);
    for i=1:nShot_this
        % get cartesian coords for this shot
        Z=zxy_0{idxparam}{i}(:,1);
        X=zxy_0{idxparam}{i}(:,2);
        Y=zxy_0{idxparam}{i}(:,3);
        
        % evaluate cylindrical coords
        R_yz=sqrt(Z.^2+Y.^2);
        H_perp=X;
        Theta=wrapToPi(atan2(Z,Y)+pi/2);  % define -Z to be 0 and measure by RH +X. range=[-pi,pi]
        
        % store cyl coords
        R0_cyl{idxparam}{i}=[R_yz,Theta,H_perp];
    end
    
    % convert cylindrical R to k-space
    k_cyl{idxparam}=R0_cyl{idxparam};
    for i=1:nShot_this
        k_cyl{idxparam}{i}(:,1)=r2k(k_cyl{idxparam}{i}(:,1));
        k_cyl{idxparam}{i}(:,3)=r2k(k_cyl{idxparam}{i}(:,3));
    end
    
end

%% crop cylindrical section
% get cylindrical sector params
cyl_dtheta=diff(configs.cylsect_theta_lims);
cyl_dtrans=2*r2k(configs.cylsect_trans_hwidth);

% indices captured in section
idx_cyl_sect=cell(configs.paramset,1);
for idxparam=1:configs.paramset
    nShot_this=size(k_cyl{idxparam},1);
    idx_cyl_sect{idxparam}=cell(nShot_this,1);
    for iShot=1:nShot_this
        idx_cyl_sect{idxparam}{iShot}=1:size(k_cyl{idxparam}{iShot},1);     % get all indices
    end
end

k_cyl_sect=cell(configs.paramset,1);
for idxparam=1:configs.paramset
    k_cyl_sect{idxparam}=cell(size(k_cyl{idxparam}));
    for i=1:nShot_this
        k_cyl_sect{idxparam}{i}=k_cyl{idxparam}{i};  % get everything
        
        % cull to angular lims
        theta_tmp=k_cyl_sect{idxparam}{i}(:,2);	% [-pi,pi]
        theta_tmp=wrapTo2Pi(theta_tmp-configs.cylsect_theta_lims(1));   % zero angle to start of lim
        ind_tmp=(theta_tmp<cyl_dtheta);   % indices lying in defined angular sector
        
        k_cyl_sect{idxparam}{i}=k_cyl_sect{idxparam}{i}(ind_tmp,:);   % cull
        idx_cyl_sect{idxparam}{i}=idx_cyl_sect{idxparam}{i}(ind_tmp);   % update captured indices
        
        % cull to transverse width
        k_trans_tmp=k_cyl_sect{idxparam}{i}(:,3);  % transverse k
        ind_tmp=(abs(k_trans_tmp)<r2k(configs.cylsect_trans_hwidth));   % indices lying in defined angular section
        
        k_cyl_sect{idxparam}{i}=k_cyl_sect{idxparam}{i}(ind_tmp,:);   % cull
        idx_cyl_sect{idxparam}{i}=idx_cyl_sect{idxparam}{i}(ind_tmp);   % update captured indices
    end
end

%% Build captured counts (cartesian)
% TODO: this is quite tedious - write function to convert cart to cyl coord systems
zxy_0_captured=cell(configs.paramset,1);    % zxy counts captured
k_zxy_captured=cell(configs.paramset,1);    % k-vects captured
for idxparam=1:configs.paramset
    nShot_this=size(zxy_0{idxparam},1);
    zxy_0_captured{idxparam}=cell(nShot_this,1);    % initialise size
    k_zxy_captured{idxparam}=cell(nShot_this,1);    % initialise size
    for iShot=1:nShot_this
        % extract captured counts
        zxy_0_captured{idxparam}{iShot}=zxy_0{idxparam}{iShot}(idx_cyl_sect{idxparam}{iShot},:);
        k_zxy_captured{idxparam}{iShot}=k_zxy{idxparam}{iShot}(idx_cyl_sect{idxparam}{iShot},:);
    end
end

% Plot - highlight captured counts from the summary point cloud
% TODO: for param modes >1
if verbose>1
    %% Real space ZXY counts
    figure(h_zxy_ff);
    hold on;
    
    plot_zxy(zxy_0_captured{1},configs.num_count_disp_more,30,'r');      % highlight all the counts captured
    
    if verbose>2
        % save plot
        fname_str='zxy_ff';
        saveas(h_zxy_ff,[configs.files.dirout,fname_str,'.png']);
        saveas(h_zxy_ff,[configs.files.dirout,fname_str,'.fig']);
    end
    
    %% k-vectors
    figure(h_k_ff);
    hold on;
    
    plot_zxy(k_zxy_captured{1},configs.num_count_disp_more,30,'r');      % highlight all the counts captured
    
    if verbose>2
        % save plot
        fname_str='k_ff';
        saveas(h_k_ff,[configs.files.dirout,fname_str,'.png']);
        saveas(h_k_ff,[configs.files.dirout,fname_str,'.fig']);
    end
end

%% Get 1D-k
k_1D_cyl=cell(configs.paramset,1);
for idxparam=1:configs.paramset
    
    k_1D_cyl{idxparam}=vertcat(k_cyl_sect{idxparam}{:});
    k_1D_cyl{idxparam}=k_1D_cyl{idxparam}(:,1);     % collate all shots
end

%% Histogram in k-space
for idxparam=1:configs.paramset 
    % get hist params
    hist_k_cyl_1D.binEdge{idxparam}=configs.hist.ed_lgk;     % get log-spaced edges
    hist_k_cyl_1D.binCent{idxparam}=sqrt(hist_k_cyl_1D.binEdge{idxparam}(1:end-1).*hist_k_cyl_1D.binEdge{idxparam}(2:end));   % GEOM avg bin centres
    
    dk_cyl_volume{idxparam}=cyl_dtrans*cyl_dtheta*(hist_k_cyl_1D.binCent{idxparam}).*diff(hist_k_cyl_1D.binEdge{idxparam});     % phase space volume in k
    
    % do log-k-spaced histogram
    hist_k_cyl_1D.N{idxparam}=histcounts(k_1D_cyl{idxparam},hist_k_cyl_1D.binEdge{idxparam});
    
    % evaluate number density
    nden_k_cyl_1D{idxparam}=(hist_k_cyl_1D.N{idxparam})./(nShot_this*detQE*dk_cyl_volume{idxparam});  % normalised for: shot, QE, phase space volume
    
    % TODO - normalise for RF outcoupling efficiency
    nden_k_cyl_1D{idxparam}=nden_k_cyl_1D{idxparam}/configs.rf_frac;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quantisation noise in histogram and number density must be resolved

% % subtract background
% nk_bgd_free=nden_k_cyl_1D{1}-nden_k_cyl_1D{2};
% nk_bgd_free(nk_bgd_free<=0)=NaN;    % handle negative density

% % clean data
% idx_isnum=~isnan(nk_bgd_free);
% nk_bgd_free_clean=nk_bgd_free(idx_isnum);
% k_clean=hist_k_cyl_1D.binCent{1}(idx_isnum);    % binCent needs to be common
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Display Histogram
if verbose>0
    h_k1D_hist=figure();
    figure(h_k1D_hist);
    hold on;
    
    for paridx=1:configs.paramset
        plot(1e-6*hist_k_cyl_1D.binCent{paridx},hist_k_cyl_1D.N{paridx},'.','MarkerSize',10); % Signal
    end
    
    set(gca,'xScale','log');
    set(gca,'yScale','log');
    axis tight;
    box on;
    grid on;
    grid minor;
    
    title('Histogram');
    xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('Number in BIN$(k)$');
%     legend({'RF sweep ON','RF sweep OFF'}); % TODO - modular for param set variability
    
%     fname_str='hist_k1D';
%     saveas(h_k1D_hist,[configs.files.dirout,fname_str,'.png']);
%     saveas(h_k1D_hist,[configs.files.dirout,fname_str,'.fig']);
end

%% Display far-field k-space density profile
if verbose>0
    h_nk_cyl_1D_log=figure();
    figure(h_nk_cyl_1D_log);
    hold on;
    
%     for idxparam=1:configs.paramset
%         plot(1e-6*hist_k_cyl_1D.binCent{idxparam},...
%             1e18*nden_k_cyl_1D{idxparam},'.','MarkerSize',10);     % scale units appropriately
%     end

    set(gca,'xScale','log');
    set(gca,'yScale','log');
    %         xlim(1e6*configs.limit.k_com);
    %         ylim(1e18*configs.limit.kdensity);
%     axis tight;
    xlim(1e-6*[min(configs.hist.ed_lgk),max(configs.hist.ed_lgk)]);
%     ylim([1e18*configs.limit.det_dark_nk,1e4]);
    ylim([1e-2,1e5]);

    grid on;
    box on;
    
%     title('1D condensate momentum profile');
    xlabel('$k$ [$\mu$m$^{-1}$]'); ylabel('$n_{\infty}(k)$ [$\mu$m$^3$]');
    
%     % plot dark count noise floor (note: must be plotted after setting axis to log scale)
%     h_darkcount=refline([0,1e18*configs.limit.det_dark_nk]);
%     h_darkcount.Color='r';
%     h_darkcount.DisplayName='Dark count';
    
%     fname_str='nk_cyl_1D_loglog';
%     saveas(h_nk_cyl_1D_log,[configs.files.dirout,fname_str,'.png']);
%     saveas(h_nk_cyl_1D_log,[configs.files.dirout,fname_str,'.fig']);
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

%% Post process core analysis data
% Smooth histogram
Nk_sm=cell(configs.paramset,1);     % smoothed hist counts
g_filt=cell(configs.paramset,1);     % gaussian filter window
for idxparam=1:configs.paramset
    % copy from raw
    Nk_sm{idxparam}=hist_k_cyl_1D.N{idxparam};
    
    %% Gaussian smoothing
    g_filt{idxparam}=gausswin(configs.smooth.nspan);    % build filter
    g_filt{idxparam}=g_filt{idxparam}/sum(g_filt{idxparam});    % normalise
    
    % apply smoothing filter
    Nk_sm{idxparam}=conv(Nk_sm{idxparam},g_filt{idxparam},'same');
    
    %% Plot
    if verbose>0
        figure(h_k1D_hist);
%         legend_str_this=sprintf('%d gauss smooth [%d]',idxparam,configs.smooth.nspan);
%         plot(1e-6*hist_k_cyl_1D.binCent{idxparam},Nk_sm{idxparam},...
%             'o','MarkerSize',10,'DisplayName',legend_str_this);
%         legend(gca,'off')
%         legend(gca,'show');
    plot(1e-6*hist_k_cyl_1D.binCent{idxparam},Nk_sm{idxparam},...
        'o','MarkerSize',10);
    end
end

% Smoothed density profile
nk_sm=cell(configs.paramset,1);     % smoothed k-density
for idxparam=1:configs.paramset
    nShot_this=size(zxy_0{idxparam},1);
    nk_sm{idxparam}=Nk_sm{idxparam}./(nShot_this*detQE*dk_cyl_volume{idxparam});  % normalised for: shot, QE, phase space volume
    % normalise for RF outcoupling efficiency
    nk_sm{idxparam}=nk_sm{idxparam}/configs.rf_frac;
    
    % Plot
    if verbose>0
        figure(h_nk_cyl_1D_log);
        legend_str_this=sprintf('%d gauss smooth [%d]',idxparam,configs.smooth.nspan);
%         plot(1e-6*hist_k_cyl_1D.binCent{idxparam},1e18*nk_sm{idxparam},...
%             '.-','MarkerSize',10,'DisplayName',legend_str_this);
        plot(1e-6*hist_k_cyl_1D.binCent{idxparam},1e18*nk_sm{idxparam},...
            '-','LineWidth',3,'Color','b');
%         legend(gca,'off')
%         legend(gca,'show');
    end
end

% get data for further analysis
k=hist_k_cyl_1D.binCent{1};     % k-bin centres (exp-spaced)
nk=nk_sm{1};            % evaluated density at k-point

% moving standard deviation
% n_sm_movingstd=3;       % TODO: this is really strange...
% nk_std=movingstd(nk,n_sm_movingstd,'central');  % smoothing standard deviation
% nk_se=nk_std/sqrt(n_sm_movingstd);

% %% Plot smoothed density profile
% if verbose>0
%     figure(h_nk_cyl_1D_log);
%     ax_nk_shaded=shadedErrorBar(1e-6*k,1e18*nk,1e18*nk_se,'-k',1);
% %     
% %     % legend
% %     ax_nk_shaded.mainLine.DisplayName='Smoothed data';
% end

%% FIT THERMAL DEPLETION density profile
%%%%%%%%%%% LOGARITHMIC
% get fitting region
I_th=zeros(1,2);        % [min,max] index lims
for idx=1:2    
    [~,I_th(idx)]=min(abs(k-configs.fit.thermal.k_lim(idx)));   % get indices to closest k-cent bin
end

% get fitting data
k_fit.thermal.k=k(I_th(1):I_th(2));
k_fit.thermal.lg_nk=log(nk(I_th(1):I_th(2)));

% call fitting routine
k_fit.thermal.fit=fitnlm(k_fit.thermal.k,k_fit.thermal.lg_nk,...
    configs.fit.thermal.fun,configs.fit.thermal.param0,...
    'CoefficientNames',configs.fit.thermal.coefname,...
    'Options',configs.fit.thermal.opt);

% Summarise fit
disp(k_fit.thermal.fit);

% build a sample of the fitted model's profile
extrap_temp=configs.fit.thermal.extrap;
k_fit.thermal.k_fit=linspace(min(k_fit.thermal.k)/extrap_temp,...
    max(k_fit.thermal.k)*extrap_temp,1000);   % indep var to evaluate fitted function
k_fit.thermal.lg_nk_fit=feval(k_fit.thermal.fit,k_fit.thermal.k_fit);  % evaluate fitted model

% Plot
figure(h_nk_cyl_1D_log); hold on;
coeffnames=k_fit.thermal.fit.CoefficientNames;
dispname=sprintf('TD: (%s,%s)=(%0.2g,%0.2g)\n SE(%0.1g,%0.1g)',...
    coeffnames{:},k_fit.thermal.fit.Coefficients.Estimate,...
    k_fit.thermal.fit.Coefficients.SE);
% plot(1e-6*k_fit.thermal.k_fit,1e18*exp(k_fit.thermal.lg_nk_fit),...
%     ':','LineWidth',2,'DisplayName',dispname);
plot(1e-6*k_fit.thermal.k_fit,1e18*exp(k_fit.thermal.lg_nk_fit),...
    ':','LineWidth',2,'Color','k');

% %%%%%%%%%%% LINEAR
% % get fitting region
% I_th=zeros(1,2);        % [min,max] index lims
% for idx=1:2    
%     [~,I_th(idx)]=min(abs(k-configs.fit.thermal.k_lim(idx)));   % get indices to closest k-cent bin
% end
% 
% % get fitting data
% k_fit.thermal.k=k(I_th(1):I_th(2));
% k_fit.thermal.nk=nk(I_th(1):I_th(2));
% 
% % call fitting routine
% k_fit.thermal.fit=fitnlm(k_fit.thermal.k,k_fit.thermal.nk,...
%     configs.fit.thermal.fun,configs.fit.thermal.param0,...
%     'CoefficientNames',configs.fit.thermal.coefname,...
%     'Options',configs.fit.thermal.opt);
% 
% % Summarise fit
% disp(k_fit.thermal.fit);
% 
% % build a sample of the fitted model's profile
% extrap_temp=configs.fit.thermal.extrap;
% k_fit.thermal.k_fit=linspace(min(k_fit.thermal.k)/ratio_extrap,...
%     max(k_fit.thermal.k)*ratio_extrap,1000);   % indep var to evaluate fitted function
% k_fit.thermal.nk_fit=feval(k_fit.thermal.fit,k_fit.thermal.k_fit);  % evaluate fitted model
% 
% % Plot
% figure(h_nk_cyl_1D_log); hold on;
% plot(1e-6*k_fit.thermal.k_fit,1e18*k_fit.thermal.nk_fit,...
%     ':','LineWidth',2,'DisplayName','Thermal depletion fit');


% update plot
% axis tight;
% legend(gca,'off')
% legend(gca,'show');

% Patch fitting region
ax_this=gca;
ylim_region=ax_this.YLim;               % Y-lim to axis lims
xlim_region=1e-6*configs.fit.thermal.k_lim;     % X-lim to fitting region

vert_region=[xlim_region flip(xlim_region);ylim_region(1) ylim_region(1) ylim_region(2) ylim_region(2)]';
face_region=1:4;

h_patch_region=patch('Faces',face_region,'Vertices',vert_region,...
    'FaceColor','red','EdgeColor','none','FaceAlpha',0.1);
uistack(h_patch_region,'bottom');   % send patch object to bottom


%% FIT QUANTUM DEPLETION density profile
% get fitting region
I_qd=zeros(1,2);        % [min,max] index lims
for idx=1:2    
    [~,I_qd(idx)]=min(abs(k-configs.fit.qd.k_lim(idx)));   % get indices to closest k-cent bin
end

% get fitting data
k_fit.QD.lg_k=log(k(I_qd(1):I_qd(2)));
k_fit.QD.lg_nk=log(nk(I_qd(1):I_qd(2)));

% call fitting routine
k_fit.QD.fit=fitnlm(k_fit.QD.lg_k,k_fit.QD.lg_nk,...
    configs.fit.qd.fun,configs.fit.qd.param0,...
    'CoefficientNames',configs.fit.qd.coefname,...
    'Options',configs.fit.qd.opt);

% Summarise fit
disp(k_fit.QD.fit);

% build a sample of the fitted model's profile
extrap_temp=configs.fit.qd.extrap;
k_fit.QD.lg_k_fit=linspace((min(k_fit.QD.lg_k))-log(extrap_temp),...
    max(k_fit.QD.lg_k)+log(extrap_temp),1000);   % indep var to evaluate fitted function
k_fit.QD.lg_nk_fit=feval(k_fit.QD.fit,k_fit.QD.lg_k_fit);  % evaluate fitted model

% Plot
figure(h_nk_cyl_1D_log); hold on;
coeffnames=k_fit.QD.fit.CoefficientNames;
dispname=sprintf('QD: (%s,%s)=(%0.2g,%0.2g)\n SE(%0.1g,%0.1g)',...
    coeffnames{:},k_fit.QD.fit.Coefficients.Estimate,...
    k_fit.QD.fit.Coefficients.SE);
% plot(1e-6*exp(k_fit.QD.lg_k_fit),1e18*exp(k_fit.QD.lg_nk_fit),...
%     '--','LineWidth',2,'DisplayName',dispname);
plot(1e-6*exp(k_fit.QD.lg_k_fit),1e18*exp(k_fit.QD.lg_nk_fit),...
    '--','LineWidth',2,'Color','k');

% update plot
% axis tight;
% legend(gca,'off')
% legend(gca,'show');

% Patch fitting region
ax_this=gca;
ylim_region=ax_this.YLim;               % Y-lim to axis lims
xlim_region=1e-6*configs.fit.qd.k_lim;     % X-lim to fitting region

vert_region=[xlim_region flip(xlim_region);ylim_region(1) ylim_region(1) ylim_region(2) ylim_region(2)]';
face_region=1:4;

h_patch_region=patch('Faces',face_region,'Vertices',vert_region,...
    'FaceColor','blue','EdgeColor','none','FaceAlpha',0.1);
uistack(h_patch_region,'bottom');   % send patch object to bottom

grid minor;
ax=gca;
ax.FontSize=13;
ax.PlotBoxAspectRatio=[1,0.85,1];
ax.TickLength=2*[0.01,0.01];


%% Save fitted k-density profile
fname_str='k_profile';

% h_nk_cyl_1D_log.Units = configs.fig.paperunits;
% h_nk_cyl_1D_log.Position = configs.fig.paperposition;
% 
% h_nk_cyl_1D_log.PaperSize = configs.fig.papersize;
% h_nk_cyl_1D_log.PaperUnits = configs.fig.paperunits;
% h_nk_cyl_1D_log.PaperPosition = configs.fig.paperposition;

% saveas(h_nk_cyl_1D_log,[configs.files.dirout,fname_str,'.png']);
% saveas(h_nk_cyl_1D_log,[configs.files.dirout,fname_str,'.fig']);
% saveas(h_nk_cyl_1D_log,[configs.files.dirout,fname_str,'_',date,'.eps'], 'psc2');     % save fig in cd
print(h_nk_cyl_1D_log,[configs.files.dirout,fname_str,'_',date],'-dsvg');


%% Report some findings
nShotCounts=shotSize(zxy_0{1})/(configs.rf_frac*configs.const.detect_qe);
N0_pred=[mean(nShotCounts),std(nShotCounts)];
fprintf('N0 =  %0.2g ± %0.2g\n',N0_pred(1),N0_pred(2));


%% Save data
% TODO - package this into a function
% TODO - if datafile already exists, MOVE it

% parse list of vars user wants to save and check against workspace if exists
varsExist=cell(size(vars_save));
varCounter=0;
for i = 1:length(vars_save)
    tVar=vars_save{i};
    if ~exist(tVar,'var')
        warning(['Variable "',tVar,'" does not exist.']);
        continue;
    end 
    varCounter=varCounter+1;
    varsExist{varCounter}=tVar;
end
varsExist=varsExist(1:varCounter);

save(configs.files.saveddata,varsExist{:});

%% END
t_main_end=toc(t_main_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
disp('===================ALL TASKS COMPLETED===================');