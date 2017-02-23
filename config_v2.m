% Config file for quantum depletion in far-field momentum distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data files iterate through predefined experimental param set
% allow for different analysis depending on param set used
%
% log-log transform of data and a linear fit to obtain power law
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% GENERAL
verbose=2;

%%% Raw data handling
% files -  data file
configs.files.path='C:\Users\HE BEC\Documents\lab\quantum-depletion\exp6\d';    % path to unindexed data file (e.g. 'a\b\datadir\$DATA_FNAME_TOKEN$')
configs.files.id=1:2696;          % file id numbers to use for analysis
configs.files.minCount=1000;     % min counts to use for analysis


% XY plane rotation to align to trap geometry
configs.rot_angle=0.61;     % angle in rad (exp default should be 0.61)

% TXY window - region of interest ( [] --> no crop )
%   XY window is applied after rotation
%   window liberally for very long k-tail
configs.window{1}=[0.50 0.64];      % T [s]
configs.window{2}=[-45e-3,40e-3];    % X [m]
configs.window{3}=[-50e-3,38e-3];    % Y [m]

% Param set
configs.paramset=1;     % TODO: number of params iterated - too specific

% BEC locator - used for accurate location of condensate
configs.bec.txy_pos=[0.5715,-3.45e-3,-10e-3];       % approx bec location (get from DLD front panel)
configs.bec.box_fwidth=[15e-3,20e-3,40e-3];     % txy full-width for box capt (be liberal-takes mean)

%%% Experimental consts 
configs.const.detect_qe=0.1;      % detector quantum efficiency
configs.const.hbar=1.05457e-34;     % hbar [m^2kg/s]
configs.const.m_He=6.646476e-27;    % mass of helium (kg)
configs.const.tof=0.416;    % TOF for free-fall from trap to DLD

% TODO - this is a duplicate - needs to be merged
configs.num_count_disp=1e5;
configs.num_count_disp_more=1e6;

%% Quantum depletion specific
%%% Angular integration over radial profile - Cylindrical sector
configs.cylsect_dtheta=deg2rad(15);
configs.cylsect_theta_lims=pi+configs.cylsect_dtheta*[-0.5,0.5];    % in-radial plane angle lims for angular averaging
configs.cylsect_trans_hwidth=k2r(0.8e6);          % transverse averaging half width [m]


%% Histogramming
configs.hist.ed_lgk=logspace(log10(0.3e6),log10(20e6),100);   % 10^X [m^-1 == 1e-6 um^-1]

% DEBUG FOR FLAT BACKGROUND
do_flat=0;

%% Smoothing
% TODO: understand this better!
configs.smooth.nspan=5;

%% Fit to large-k tail
% Polylog function for non-integer order: polylog2

%%% Thermal depletion
% fitting region
configs.fit.k_lim_thermal=[2.5e6 6e6];      % bound k to use for fit [m^-1]
% fitting function
configs.fit.fun_thermal='y~Nth*polylog2(3/2,exp(-6.0597e-20*(x1^2)/Ta))/(4.4870e-29*(Ta^(3/2)))';    % negative-power function
configs.fit.fun_coefname_thermal={'Nth','Ta'};         % function coefficient names
% initial conditions
configs.fit.param0_thermal=[1e3,1e-8];   % (A, alpha)
% fit options
configs.fit.opt_thermal=statset('TolFun',1e-15,...
    'TolX',1e-15,...
    'MaxIter',1e6,...
    'UseParallel',1,...
    'Display','off');

%%% Quantum depletion - loglog transformed linear
% fitting region
configs.fit.k_lim=[6e6 15e6];      % k-region of data to fit

% fitting function and init conditions
% configs.fit.fun_negpowk='y~A-alpha*x1';     % linearised fitting function
% configs.fit.fun_coefname={'A','alpha'};    	% function coefficient names
% configs.fit.param0=[10,4.0];   % (A, alpha)

configs.fit.fun_negpowk='y~log(C_inf/248.0502)-alpha*x1';     % linearised fitting function
configs.fit.fun_coefname={'C_inf','alpha'};    	% function coefficient names
configs.fit.param0=[(2*pi)^3*exp(10),4.0];  % (C_inf, alpha)

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

configs.limit.det_dark_nk=6.5e-3*1e-18;   % dark ct density in k [m^3]
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