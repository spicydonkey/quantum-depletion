% k-space Thomas-Fermi distribution
% params: [n0, c]

function lg_n_tf = log_tf_dist(params, k)
    n0=params(1);
    c=params(2);
%     
%     m_he=6.6464e-27;    % mass of helium [kg]
%     h=6.6261e-34;      	% Planck's constant [m^2.kg/s]
%     k_B=1.3807e-23;   	% Boltzmann constant [m^2.kg.s^-2.K^-1]
%         
    n_tf=n0-c*(k.^2);
    n_tf(n_tf<0)=0;     % negative density to zero
    
    lg_n_tf=log(n_tf);
    
end