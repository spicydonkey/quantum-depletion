% TXY_FLAT is from a flat background (3D linearly spaced points in a similar
% range as experiment)
% NOTE: TXY_FLAT is origin centred

% Lims
LIM{1}=[-0.1,0.1];          % T[s]
LIM{2}=[-40e-3,40e-3];      % X[m]
LIM{3}=[-40e-3,40e-3];      % Y[m]

if ~exist('N_PER_SHOT','var')
    N_PER_SHOT=10;
end
if ~exist('N_SHOT','var')
    N_SHOT=1000;
end

TXY_FLAT=cell(N_SHOT,1);

for iShot=1:N_SHOT
    % monte-carlo
    TXY_FLAT{iShot}=rand(N_PER_SHOT,3);
    
    for idx=1:3
        TXY_FLAT{iShot}(:,idx)=TXY_FLAT{iShot}(:,idx)*diff(LIM{idx})+LIM{idx}(1);
    end
end

% Plot
figure();
plot_zxy(TXY_FLAT);
view(3);
title('Simulated flat background');
xlabel('X [m]');ylabel('Y [m]');zlabel('T [s]');