%% Configs
configs.files.path='\\AMPLPC29\He BEC Archive\EXPERIMENT-DATA\QuantumDepletion\rf_oc_frac\vpp_300\d';

configs.load.id=1:7;
configs.load.mincount=0;            % min counts in window - 0 for no min
configs.load.maxcount=Inf;          % max counts in window - Inf for no max
configs.load.rot_angle=0.61;

configs.load.window{1}=[0.56,0.58];      % T [s]
configs.load.window{2}=[-35e-3,35e-3];    % X [m]
configs.load.window{3}=[-35e-3,35e-3];    % Y [m]

% define region for outcoupled atom (mf=0)
vz=9.81*0.416;
boxLim={[0.56,0.58]*vz,[-0.01,0.005],[-35e-3,35e-3]};