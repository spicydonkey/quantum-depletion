% Config file for quantum depletion in far-field momentum distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data files iterate through predefined experimental param set
% allow for different analysis depending on param set used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% GENERAL
verbose=2;

%%% Raw data handling
% files -  data file
configs.files.path='C:\Users\HE BEC\Documents\lab\quantum-depletion\exp4\d';    % path to unindexed data file (e.g. 'a\b\datadir\$DATA_FNAME_TOKEN$')
configs.files.id=1:1500;          % file id numbers to use for analysis
configs.files.minCount=100;     % min counts to use for analysis


% XY plane rotation to align to trap geometry
configs.rot_angle=0.61;     % angle in rad (exp default should be 0.61)

% TXY window - region of interest ( [] --> no crop )
%   XY window is applied after rotation
%   window liberally for very long k-tail
configs.window{1}=[0.50 0.64];      % T [s]
configs.window{2}=[-45e-3,40e-3];    % X [m]
configs.window{3}=[-50e-3,38e-3];    % Y [m]


% Param set
configs.paramset=2;     % TODO: number of params iterated


% BEC locator - used for accurate location of condensate
configs.bec.txy_pos=[0.5715,-3.45e-3,-10e-3];       % approx bec location (get from DLD front panel)
configs.bec.box_fwidth=[15e-3,20e-3,40e-3];     % txy full-width for box capt (be liberal-takes mean)


%%% Experimental consts 
configs.const.detect_qe=0.1;      % detector quantum efficiency
configs.const.hbar=1.05457e-34;     % hbar [m^2kg/s]
configs.const.m_He=6.646476e-27;    % mass of helium (kg)
configs.const.tof=0.416;    % TOF for free-fall from trap to DLD


%% Quantum depletion specific
%%% Angular integration over radial profile - Cylindrical sector
configs.cylsect_theta_lims=pi+[-pi/8,pi/8];    % in-radial plane angle lims for angular averaging
configs.cylsect_trans_hwidth=4e-3;          % transverse averaging half width [m]

%% 1D slice
% TODO - document usage properly
% counts captured along a 1D line-slice from well "below" condensate to centre
configs.slice.mincount=500;     % minimum count in 1D slice to pass
configs.slice.cyl_orient=1;     % slice through Z-axis (1:Z,2:X,3:Y)
configs.slice.cyl_rad=1e-3;     % cyl radius [m]
configs.slice.cyl_hgt=280e-3;    % cyl height [m]

% build cylinder dim
configs.slice.cyl_dim=[configs.slice.cyl_rad,configs.slice.cyl_hgt];

% build cylinder centre such that slice ends at condensate centre
configs.slice.cyl_cent=zeros(1,3);
configs.slice.cyl_cent(configs.slice.cyl_orient)=-0.5*configs.slice.cyl_hgt;

%%% Angular averaging
% Rotate around X-direction (radial plane is YZ)
% TODO - feature isn't implemented - just uses {1}
configs.axial_rot_angle{1}=linspace(-pi/8,pi/8,20);     % point -Z
configs.axial_rot_angle{2}=pi+linspace(-pi/8,pi/8,20);  % point +Z


%% Background density
% configs.do_bgd_calc=1;
% +Z direction
configs.bgd_cyl_dir{1}=1;
configs.bgd_cyl_orient{1}=1;       
configs.bgd_disp{1}=[0,0.027,0.027];    % T,X,Y disp to do background analysis
configs.bgd_cyl_rad{1}=5e-3;
configs.bgd_cyl_hgt{1}=k2r(20*1e6);     % cylinder oriented in +Z direction

% -Z direction
configs.bgd_cyl_dir{2}=-1;     % cylinder oriented in -Z direction
configs.bgd_cyl_orient{2}=1;       % cylinder oriented in -Z direction
configs.bgd_disp{2}=[0,0.027,0.027];    % T,X,Y disp to do background analysis
configs.bgd_cyl_rad{2}=5e-3;
configs.bgd_cyl_hgt{2}=k2r(20*1e6);    

% build cylinder params
for i=1:length(configs.bgd_cyl_orient)
    configs.bgd_cyl_dim{i}=[configs.bgd_cyl_rad{i},configs.bgd_cyl_hgt{i}];
    configs.bgd_cyl_cent{i}=configs.bgd_disp{i};
    configs.bgd_cyl_cent{i}(configs.bgd_cyl_orient{i})=configs.bgd_cyl_cent{i}(configs.bgd_cyl_orient{i})+0.5*configs.bgd_cyl_dir{i}*configs.bgd_cyl_hgt{i};
end


%% Histogramming
configs.hist.nbin=100;   % used for real space and linear-k distribution
% TODO - (log-spaced) bin edges to use in k-space --> simplifies anisotropic
% analysis
configs.hist.ed_lgk=logspace(5,7.3,1000);   % 10^X [m^-1 == 1e-6 um^-1]


%% Fit to large-k tail
%% TODO
% %%% Thermal depletion
% % fitting region
% configs.fit.k_lim_thermal=[2.5e6 6e6];      % bound k to use for fit [m^-1]
% % fitting function
% configs.fit.fun_negpowk_thermal='y~(x1^())';    % negative-power function
% configs.fit.fun_coefname_thermal={'A','alpha'};         % function coefficient names
% % initial conditions
% configs.fit.param0_thermal=[1e10,4];   % (A, alpha)
% % fit options
% configs.fit.opt_thermal=statset('TolFun',1e-15,...
%     'TolX',1e-15,...
%     'MaxIter',1e6,...
%     'UseParallel',1,...
%     'Display','off');

%%% Quantum depletion
% fitting region
configs.fit.k_min=8e6;      % lower bound k to use for fit [m^-1]
% fitting function
configs.fit.fun_negpowk='y~A*(x1^(-alpha))';    % negative-power function
configs.fit.fun_coefname={'A','alpha'};         % function coefficient names
% initial conditions
configs.fit.param0=[1e10,4.0];   % (A, alpha)
% fit options
configs.fit.opt=statset('TolFun',1e-30,...
    'TolX',1e-30,...
    'MaxIter',1e8,...
    'UseParallel',1,...
    'Display','off');


%% Limits
%%% Scale
% size of cloud limited by detector resolution (0.1mm) and BEC (68mm=10um-1)
configs.limit.r_com=1e-3*[0.1,68];      % [m]
configs.limit.k_com=1e-6*[0.1,10];      % [m^-1]

%%% Density 
% TODO set more accurately - currently rounded to power of 10
% evaluated from dark counts (~6.5e-3 um3) and saturation (~100 mm-^3)
configs.limit.det_dark_nr=1e9*1e-5;  	% dark ct density in r [m^-3]
configs.limit.det_sat_nr=1e9*1e3;       % sat ct density in r

configs.limit.det_dark_nk=1e-18*1e-3;   % dark ct density in k [m^3]
configs.limit.det_sat_nk=1e-18*1e5;     % sat ct density in k

% create limits for plotting
configs.limit.rdensity=[configs.limit.det_dark_nr,configs.limit.det_sat_nk];
configs.limit.kdensity=[configs.limit.det_dark_nk,configs.limit.det_sat_nk];


%% ALGORITHM CONFIGS
% DO NOT ADJUST
configs.files.saveddata=[configs.files.path,'_data.mat'];     % path to store saved data
configs.files.archive=[configs.files.path,'_archive'];   % dir to archive folder
configs.files.dirout=[configs.files.path,'_output\'];      % output directory

% create dir if necessary
if ~exist(configs.files.archive,'dir')
    warning('Archive directory %s does not exist. Creating directory...',configs.files.archive);
    mkdir(configs.files.archive);
end
if ~exist(configs.files.dirout,'dir')
    warning('Output directory %s does not exist. Creating directory...',configs.files.dirout);
    mkdir(configs.files.dirout);
end