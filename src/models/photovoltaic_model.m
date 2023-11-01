function [i_ps, v_ps, p_ps] = photovoltaic_model(irradiation)
    % Photovoltaic panel model
    
    % Declare constants
    Q = 1.6e-19; % Elementary charge
    K = 1.38e-23; % Boltzmann constant
    
    % Declare photovoltaic parameters
    N_ps = 8;            % Number of panels in parallel
    N_ss = 300;          % Number of panels in series
    T_ps = 298;          % Temperature
    Tr = 298;            % Reference temperature
    Isc = 3.27;          % Short circuit current at Tr
    Kl = 0.0017;         % Short circuit current temperature coeff
    Ior = 2.0793e-6;     % Ior - Irs at Tr
    Ego = 1.1;           % Band gap energy of the semiconductor
    A = 1.6;             % Factor. cell deviation from de ideal pn junction   

    % Intermediate photovoltaic variables
    Vt = K*T_ps/Q;
    Iph = (Isc+Kl*(T_ps-Tr))*irradiation;
    Irs = Ior*(T_ps/Tr)^3*exp(Q*Ego*(1/Tr-1/T_ps)/(K*A));

    % Algebraic equations
    v_ps = (N_ss*Vt*A*(lambertw_approx(exp(1)*(Iph/Irs+1))-1));
    i_ps = N_ps*(Iph-Irs*(exp(v_ps/(N_ss*Vt*A))-1));
    p_ps = v_ps*i_ps;
end