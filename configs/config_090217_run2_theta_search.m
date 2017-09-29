% Config file for quantum depletion in far-field momentum distribution
% mf=0 atoms

%%% GENERAL
verbose=1;

%%% Raw data handling
% files -  data file
configs.files.path='C:\Users\HE BEC\Documents\lab\quantum-depletion\exp2\d';    % path to unindexed data file (e.g. 'a\b\datadir\$DATA_FNAME_TOKEN$')
configs.files.id=1:2500;          % file id numbers to use for analysis
configs.files.minCount=100;     % min counts to use for analysis


% XY plane rotation to align to trap geometry
configs.rot_angle=0.61;     % angle in rad (exp default should be 0.61)

% TXY window - region of interest ( [] --> no crop )
%   XY window is applied after rotation
%   window liberally for very long k-tail
configs.window{1}=[0.54 0.595];      % T [s]
configs.window{2}=[-45e-3,40e-3];    % X [m]
configs.window{3}=[-50e-3,38e-3];    % Y [m]


% BEC locator - used for accurate location of condensate
configs.bec.txy_pos=[0.567,-2e-3,-10e-3];       % approx bec location (get from DLD front panel)
configs.bec.box_fwidth=[15e-3,20e-3,40e-3];     % txy full-width for box capt (be liberal-takes mean)


%%% Experimental consts 
configs.const.detect_qe=0.1;      % detector quantum efficiency
configs.const.hbar=1.05457e-34;     % hbar [m^2kg/s]
configs.const.m_He=6.646476e-27;    % mass of helium (kg)
configs.const.tof=0.430;    % TOF for free-fall from trap to DLD


%% Quantum depletion specific
%%% Angular integration over radial profile - Cylindrical sector
configs.cylsect_theta_lims=[-pi/4,pi/4];    % in-radial plane angle lims for angular averaging
configs.cylsect_trans_hwidth=5e-3;          % transverse averaging half width [m]

%%% 1D slice
% counts captured along a 1D line-slice from well "below" condensate to centre
configs.slice.cyl_orient=1;     % slice through Z-axis (1:Z,2:X,3:Y)
configs.slice.cyl_rad=1e-3;     % cyl radius [m]
configs.slice.cyl_hgt=80e-3;    % cyl height [m]
configs.slice.mincount=500;     % minimum count in 1D slice to pass

% build cylinder dim
configs.slice.cyl_dim=[configs.slice.cyl_rad,configs.slice.cyl_hgt];

% build cylinder centre such that slice ends at condensate centre
configs.slice.cyl_cent=zeros(1,3);
configs.slice.cyl_cent(configs.slice.cyl_orient)=-0.5*configs.slice.cyl_hgt;

%%% Histogramming
configs.hist.nbin=100;   % used for real space and linear-k distribution
% TODO - (log-spaced) bin edges to use in k-space --> simplifies anisotropic
% analysis
configs.hist.ed_lgk=logspace(5,7,300);   % 10^X [m^-1 == 1e-6 um^-1]

%%% Angular averaging
% Rotate around X-direction (radial plane is YZ)
configs.axial_rot_angle=linspace(-pi,pi,200);     % angles to perform 1D analysis

%%% Fit to large-k tail
% fitting region
configs.fit.k_min=3.5e6;      % lower bound k to use for fit [m^-1]
% fitting function
configs.fit.fun_negpowk='y~A*(x1^(-alpha))';    % negative-power function
configs.fit.fun_coefname={'A','alpha'};         % function coefficient names
% initial conditions
configs.fit.param0=[1e10,4];   % (A, alpha)
% fit options
configs.fit.opt=statset('TolFun',1e-15,...
    'TolX',1e-15,...
    'MaxIter',1e6,...
    'UseParallel',1,...
    'Display','off');


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