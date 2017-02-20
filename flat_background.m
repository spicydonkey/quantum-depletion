% TXY_FLAT is from a flat background (3D linearly spaced points in a similar
% range as experiment)

% Lims
LIM{1}=[-0.075,0.075];   % T[s]
LIM{2}=[-40e-3,40e-3];   % X[m]
LIM{3}=[-40e-3,40e-3];   % Y[m]

N_PER_SHOT=10;
N_SHOT=3000;

TXY_FLAT=cell(N_SHOT,1);

for iShot=1:N_SHOT
    % monte-carlo
    TXY_FLAT{iShot}=rand(N_PER_SHOT,3);
    
    for idx=1:3
        TXY_FLAT{iShot}(:,idx)=TXY_FLAT{iShot}(:,idx)*diff(LIM{idx})+LIM{idx}(1);
    end
end

