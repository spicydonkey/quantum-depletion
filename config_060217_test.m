% Config file for quantum depletion in far-field momentum distribution
% mf=0 atoms

%%% GENERAL
verbose=2;

%%% Raw data handling
% files -  data file
configs.files.path='\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\QuantumDepletion\Run3\d';    % path to unindexed data file (e.g. 'a\b\datadir\$DATA_FNAME_TOKEN$')
configs.files.id=1:3500;          % file id numbers to use for analysis
configs.files.minCount=100;     % min counts to use for analysis


% XY plane rotation to align to trap geometry
configs.rot_angle=0.61;     % angle in rad (exp default should be 0.61)

% TXY window - region of interest ( [] --> no crop )
%   XY window is applied after rotation
%   window liberally for very long k-tail
configs.window{1}=[0.54 0.585];      % T [s]
configs.window{2}=[-35e-3,30e-3];    % X [m]
configs.window{3}=[-40e-3,28e-3];    % Y [m]


%%% Experimental consts 
configs.const.detect_qe=0.1;      % detector quantum efficiency
configs.const.hbar=1.05457e-34;     % hbar [m^2kg/s]
configs.const.m_He=6.646476e-27;    % mass of helium (kg)
configs.const.tof=0.430;    % TOF for free-fall from trap to DLD


%% Quantum depletion specific
%%% 1D slice
% counts captured along a 1D line-slice from well "below" condensate to centre
configs.slice.cyl_orient=1;     % slice through Z-axis (1:Z,2:X,3:Y)
configs.slice.cyl_rad=5e-3;     % cyl radius [m]
configs.slice.cyl_hgt=70e-3;    % cyl height [m]

% build cylinder dim
configs.slice.cyl_dim=[configs.slice.cyl_rad,configs.slice.cyl_hgt];

% build cylinder centre such that slice ends at condensate centre
configs.slice.cyl_cent=zeros(1,3);
configs.slice.cyl_cent(configs.slice.cyl_orient)=-0.5*configs.slice.cyl_hgt;

%%% HISTOGRAMMING
configs.nbin=100;


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