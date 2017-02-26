% k-space Bose distribution
% params: [N_THERMAL, T_APPARENT]

function lg_n_bose = log_bose_dist(params, k)
    Nth=params(1);
    Ta=params(2);
%     
%     m_he=6.6464e-27;    % mass of helium [kg]
%     h=6.6261e-34;      	% Planck's constant [m^2.kg/s]
%     k_B=1.3807e-23;   	% Boltzmann constant [m^2.kg.s^-2.K^-1]
%     
%     EQNCONST=h/sqrt(2*pi*m_he*k_B);
%     
%     L_dB=EQNCONST/sqrt(Ta);     % de-Broglie length
%     
%     lg_n_bose=log(Nth*g_bose(-(k.^2)*L_dB^2/(4*pi)) / (1.202*(2*pi/L_dB)^3));
    
    lg_n_bose=log(Nth*g_bose(exp(-6.0597e-20*(k.^2)/Ta))/(4.4870e29*(Ta^(3/2))));
end