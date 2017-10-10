% Characterise RF sweep outcoupling fraction

%% Load config
path_config='C:\Users\HE BEC\Documents\MATLAB\quantum-depletion\configs\config_rf_0v9';
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

eff_rf_outcoupling=numCountsRF./numCountsTot;
fprintf('RF outcoupling efficiency = %0.3g (%0.1g)\n',mean(eff_rf_outcoupling),std(eff_rf_outcoupling));

% plot
figure();
plot_zxy(zxy,[],1,'k');
hold on;
plot_zxy(zxy_mf0,[],5,'r');

xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');