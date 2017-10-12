%% Configs
configs.files.path='\\AMPLPC29\He BEC Archive\EXPERIMENT-DATA\QuantumDepletion\rf_oc_frac\vpp_1300_zpush\d';

configs.load.id=1:60;
configs.load.mincount=100;            % min counts in window - 0 for no min
configs.load.maxcount=Inf;          % max counts in window - Inf for no max
configs.load.rot_angle=0.61;

configs.load.window{1}=[0.48,0.56];      % T [s]
configs.load.window{2}=[-35e-3,35e-3];    % X [m]
configs.load.window{3}=[-35e-3,35e-3];    % Y [m]

% define region for outcoupled atom (mf=0)
vz=9.81*0.416;

% box regions for triplets
% IDX   |   MF
%-------+----------
% 1     |   1
% 2     |   0
% 3     |   -1
mf=[1,0,-1];
boxLim{1}={[0.48,0.51]*vz,[-35e-3,35e-3],[-35e-3,35e-3]};  
boxLim{2}={[0.51,0.528]*vz,[-35e-3,35e-3],[-35e-3,35e-3]};
boxLim{3}={[0.528,0.56]*vz,[-35e-3,35e-3],[-35e-3,35e-3]};