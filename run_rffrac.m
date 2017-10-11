% Characterise RF sweep outcoupling fraction

%% Load config
path_config='C:\Users\HE BEC\Documents\MATLAB\quantum-depletion\configs\config_rf_zpush';
run(path_config);

%% Load data
[txy,fout]=load_txy(configs.files.path,configs.load.id,...
    configs.load.window,configs.load.mincount,configs.load.maxcount,...
    configs.load.rot_angle,0,1,0);

zxy=txy2zxy(txy);
clearvars txy;


%% evaluate RF outcoupling fraction
zxy_mf0=cellfun(@(x) boxcull(x,boxLim),zxy,'UniformOutput',false);

numCountsTot=shotSize(zxy);       
numCountsRF=shotSize(zxy_mf0);    

[numCountsTotSorted,Isort]=sort(numCountsTot);
numCountsRFSorted=numCountsRF(Isort);

eff_rf_outcoupling=numCountsRF./numCountsTot;
eff_rf_outcoupling_sorted=numCountsRFSorted./numCountsTotSorted;
fprintf('RF outcoupling efficiency = %0.3g (%0.1g)\n',mean(eff_rf_outcoupling),std(eff_rf_outcoupling));


%% plot all counts
figure();
plot_zxy(zxy,[],1,'k');
hold on;
plot_zxy(zxy_mf0,[],5,'r');

xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');

%% plot efficiency vs total number
figure();
plot(numCountsTotSorted,eff_rf_outcoupling_sorted,'ko','MarkerSize',6);

box on;
ax=gca;
xlabel('Total number detected');
ylabel('$\eta_{RF}$');
ax.FontSize=12;