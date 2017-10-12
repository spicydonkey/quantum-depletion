% Characterise RF sweep outcoupling fraction

%% Load config
path_config='C:\Users\HE BEC\Documents\MATLAB\quantum-depletion\configs\config_rf_zpush';
run(path_config);

%% Load data
[txy,fout]=load_txy(configs.files.path,configs.load.id,...
    configs.load.window,configs.load.mincount,configs.load.maxcount,...
    configs.load.rot_angle,0,1,0);

zxy=txy2zxy(txy);


%% evaluate RF outcoupling fraction
numCountsTot=shotSize(zxy); 

% evaluate for each MF triplet state
zxy_mf=cell(1,3);
numCountsMf=cell(1,3);
eff_mf=cell(1,3);
for ii=1:3
    zxy_mf{ii}=cellfun(@(x) boxcull(x,boxLim{ii}),zxy,'UniformOutput',false);
    numCountsMf{ii}=shotSize(zxy_mf{ii});   
    eff_mf{ii}=numCountsMf{ii}./numCountsTot;
end
idx_mf0=2;
fprintf('RF outcoupling efficiency = %0.3g (%0.1g)\n',mean(eff_mf{idx_mf0}),std(eff_mf{idx_mf0}));


%% evaluate detector count rate
t_lims=configs.load.window{1};
dt_bin=0.1e-4;
% n_bin_tof=100;
n_bin_tof=round(diff(t_lims)/dt_bin);

% config gaussian filter
g_hsize=round(n_bin_tof/100);
g_sigma=round(g_hsize/5);
g_filt=gaussFilter(g_hsize,g_sigma);

t_ed=linspace(t_lims(1),t_lims(2),n_bin_tof);
t_ct=t_ed(1:end-1)+0.5*diff(t_ed);

countRate=cellfun(@(x) histcounts(x(:,1),t_ed),txy,'UniformOutput',false);
countRate=cellfun(@(x) x/dt_bin,countRate,'UniformOutput',false);
countRateFilt=cellfun(@(x) conv(x,g_filt,'same'),countRate,'UniformOutput',false);
maxCountRate=cellfun(@(x) max(x),countRateFilt);

% example TOF plot
n_shot_plot=5;
shotsToPlot=randperm(numel(txy),n_shot_plot);

cc=distinguishable_colors(n_shot_plot);

figure;
for ii=1:n_shot_plot
    thisShotId=shotsToPlot(ii);
    hold on;
    plot(t_ct,countRateFilt{thisShotId}*1e-6,...
        'Color',cc(ii,:),'LineWidth',1,...
        'DisplayName',num2str(numCountsTot(thisShotId)));
    hold off;
end


ax=gca;
ax.FontSize=12;
box on;
xlabel('TOF [s]');
ylabel('Detector count rate [MHz]');
leg=legend('show');
leg.Title.String='Total number';

%% plot all counts
figure();

cc='kbr';
for ii=1:3
    hold on;
    plot_zxy(zxy_mf{ii},[],5,cc(ii));
    hold off;
end

box on;
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');

%% plot efficiency vs total number
figure();

cc=viridis(4);
pMarker={'^','o','s'};
for ii=1:3
    hold on;
    plot(numCountsTot*1e-3,eff_mf{ii},...
        'Color',cc(ii,:),'LineStyle','none','Marker',pMarker{ii},'MarkerSize',6,...
        'DisplayName',num2str(mf(ii)));
    hold off;
end

box on;
ax=gca;
xlabel('Total number detected ($10^3$)');
ylabel('Number fraction');
ax.FontSize=12;
leg=legend('show');
leg.Title.String='$m_F$';

%% plot efficiency vs max count rate
figure();

cc=viridis(4);
pMarker={'^','o','s'};
for ii=1:3
    hold on;
    plot(maxCountRate*1e-6,eff_mf{ii},...
        'Color',cc(ii,:),'LineStyle','none','Marker',pMarker{ii},'MarkerSize',6,...
        'DisplayName',num2str(mf(ii)));
    hold off;
end

box on;
ax=gca;
xlabel('Maximum count rate [MHz]');
ylabel('Number fraction');
ax.FontSize=12;
leg=legend('show');
leg.Title.String='$m_F$';