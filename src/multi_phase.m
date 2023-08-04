%% Multi-phase 
import casadi.*

clearvars
load("inputs\data.mat")
addpath('models')

%% Problem constants

Tf = 1440;   % Final time (min)
N = 14;      % Number of control intervals per phase

M_0 = 164;     % Initial volume of hydrogen (g)
M_min = 49.2;  % Minimum volume of hydrogen (g)
M_max = 205;   % Maximum volume of hydrogen (g)
I_e_0 = 30;    % Initial current (A)
I_e_min = 30;  % Minimum current (A)
I_e_max = 100; % Maximum current (A)
I_e_N = 10;    % Stadby current (A)
tp1_0 = 480;   % Initial final time phase 1 (min)
tp2_0 = 960;   % Initial final time phase 2 (min)

%% Problblem variables

% Declare variables
v_h2 = MX.sym('v_h2'); % State - Volume of hydrogen
i_el = MX.sym('i_el'); % Control - Electrical current in electrolyzer
time = MX.sym('time'); % Time

%% Interpolate input data

Irradiation = casadi.interpolant('Irradiation', 'bspline', ...
    {t_file}, irradiation');
HydrogenDemand = casadi.interpolant('HydrogenDemand', 'bspline', ...
    {t_file}, fdemanda');

%% Models

% Dynamic
[f_h2, v_el, p_el] = electrolyzer_model(i_el);
[i_ps, v_ps, p_ps] = photovoltaic_model(Irradiation(time));
dv_h2_on = f_h2 - HydrogenDemand(time);
dv_h2_st = -HydrogenDemand(time);

% Lagrange cost function
l_on = (p_el - p_ps)^2;
l_st = (p_el - p_ps)^2;

% Right-hand side
f_on = casadi.Function('f', ...
    {v_h2, i_el, time}, {dv_h2_on, l_on}, ...
    {'x', 'u', 't'}, {'x_dot', 'L'});
f_st = casadi.Function('f', ...
    {v_h2, i_el, time}, {dv_h2_st, l_st}, ...
    {'x', 'u', 't'}, {'x_dot', 'L'});

%% Create the NPL

% Time grid
nP = 3; % number of phases
nGrid = 2*N-1; % number of points per phase
tau = linspace(0, 1, N);

X = casadi.MX.sym('X',nP*nGrid);
U = casadi.MX.sym('U',nP*nGrid);
T = [casadi.MX(0); casadi.MX.sym('T',nP-1); casadi.MX(Tf)];

L = 0;
g = MX(nP*(nGrid-1) + 3, 1);

% phase 1 - standby
X1 = X(1:nGrid);
U1 = U(1:nGrid);
for k=1:2:nGrid-1
    i = ceil(k/2);
    
    ti = T(1) + (T(2) - T(1))*tau(i);
    tf = T(1) + (T(2) - T(1))*tau(i+1);
    hk = (T(2) - T(1))*(tau(i+1) - tau(i));
    
    [f_k_0, w_k_0] = f_st(X(k), U(k), ti);
    [f_k_1, w_k_1] = f_st(X(k+1), U(k+1), ti + hk/2);
    [f_k_2, w_k_2] = f_st(X(k+2), U(k+2), tf);
    
    L = L + hk*(w_k_0 + 4*w_k_1 + w_k_2)/6;
    
    g(k) =  X(k+2) - X(k) - hk*(f_k_0 + 4*f_k_1 + f_k_2)/6;
    g(k+1) =  X(k+1) - (X(k+2) + X(k))/2 - hk*(f_k_0 - f_k_2)/8;
end

% phase 2 - active
X2 = X(nGrid+1:2*nGrid);
U2 = U(nGrid+1:2*nGrid);
for k=1:2:nGrid-1
    i = ceil(k/2);
    
    ti = T(2) + (T(3) - T(2))*tau(i);
    tf = T(2) + (T(3) - T(2))*tau(i+1);
    hk = (T(3) - T(2))*(tau(i+1) - tau(i));
    
    [f_k_0, w_k_0] = f_on(X2(k), U2(k), ti);
    [f_k_1, w_k_1] = f_on(X2(k+1), U2(k+1), ti + hk/2);
    [f_k_2, w_k_2] = f_on(X2(k+2), U2(k+2), tf);
    
    L = L + hk*(w_k_0 + 4*w_k_1 + w_k_2)/6;
    
    j = k + nGrid-1;
    g(j) =  X2(k+2) - X2(k) - hk*(f_k_0 + 4*f_k_1 + f_k_2)/6;
    g(j+1) =  X2(k+1) - (X2(k+2) + X2(k))/2 - hk*(f_k_0 - f_k_2)/8;
end

% phase 3 - standby
X3 = X(2*nGrid+1:3*nGrid);
U3 = U(2*nGrid+1:3*nGrid);
for k=1:2:nGrid-1
    i = ceil(k/2);
    
    ti = T(3) + (T(4) - T(3))*tau(i);
    tf = T(3) + (T(4) - T(3))*tau(i+1);
    hk = (T(4) - T(3))*(tau(i+1) - tau(i));
    
    [f_k_0, w_k_0] = f_st(X3(k), U3(k), ti);
    [f_k_1, w_k_1] = f_st(X3(k+1), U3(k+1), ti + hk/2);
    [f_k_2, w_k_2] = f_st(X3(k+2), U3(k+2), tf);
    
    L = L + hk*(w_k_0 + 4*w_k_1 + w_k_2)/6;
    
    j = k + 2*(nGrid-1);
    g(j) =  X3(k+2) - X3(k) - hk*(f_k_0 + 4*f_k_1 + f_k_2)/6;
    g(j+1) =  X3(k+1) - (X3(k+2) + X3(k))/2 - hk*(f_k_0 - f_k_2)/8;
end

g(nP*(nGrid-1) + 1) = X1(end) - X2(1);
g(nP*(nGrid-1) + 2) = X2(end) - X3(1);
g(nP*(nGrid-1) + 3) = T(3) - T(2);

lbg = zeros(length(g),1);
ubg = zeros(length(g),1);
lbg(end) = 0; 
ubg(end) = Tf; 

%% Set Path constrains and initial conditions

x = [X;U;T(2:3)];

lbx = [M_min*ones(nP*nGrid,1); ...
    I_e_N*ones(nGrid,1); ...
    I_e_min*ones(nGrid,1); ...
    I_e_N*ones(nGrid,1); ...
    0;0];
ubx = [M_max*ones(nP*nGrid,1); ...
    I_e_N*ones(nGrid,1); ...
    I_e_max*ones(nGrid,1); ...
    I_e_N*ones(nGrid,1); ...
    Tf; Tf];

% Set the initial condition for the state
lbx(1) = M_0;
ubx(1) = M_0;

% Set guess
xGuess = (M_max - M_min)/2;
uGuess = (I_e_max - I_e_min)/2;
x0 = [xGuess*ones(length(X),1); ...
    I_e_N*ones(length(U1),1); ...
    uGuess*ones(length(U2),1); ...
    I_e_N*ones(length(U3),1); ...
    tp1_0; tp2_0];

%% Solve

nlp = struct('x', x, 'f', L, 'g', g);
solver = nlpsol('solver', 'ipopt', nlp);

solution = solver('x0', x0, 'lbx', lbx, 'ubx', ubx,...
    'lbg', lbg, 'ubg', ubg);

%% Get solution

w_opt = (full(solution.x));
x_opt = w_opt(1:(end-2)/2);
u_opt = w_opt((end-2)/2+1:(end-2));
t_p = w_opt(end-1:end);
t_opt = [linspace(0,t_p(1), nGrid), ...
    linspace(t_p(1),t_p(2), nGrid), ...
    linspace(t_p(2),1440, nGrid)];

fprintf('\nCost: %f W^2\n',full(solution.f))
fprintf('Grid energy: %4.2f W\n',sqrt(full(solution.f)))

%% Simulate solution

nOpt = length(u_opt);
f_h2_prd_sim =  zeros(nOpt,1);
v_el_sim =  zeros(nOpt,1);
p_el_sim =  zeros(nOpt,1);
p_ps_sim = zeros(nOpt,1);
for i=1:nOpt
    [f_h2_prd_sim(i), v_el_sim(i), p_el_sim(i)] = electrolyzer_model(u_opt(i));
    [~, ~, p_ps_sim(i)] = photovoltaic_model(full(Irradiation(t_opt(i))));
end

%% Plots

figure()
subplot(3,1,1);
plot(t_opt/60, p_el_sim, 'LineWidth', 1.5, 'Color', '#0152a1')
hold on
grid on
plot(t_opt/60, p_ps_sim, 'LineWidth', 1.5, 'Color', '#ff8700')
legend('P_E','P_S')
ylabel('Power [W]')
xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])

subplot(3,1,2);
plot(t_opt/60, f_h2_prd_sim, 'LineWidth', 1.5, 'Color', '#0152a1')
hold on
grid on
plot((0:1:1440)/60,full(HydrogenDemand(0:1:1440)), 'LineWidth', 1.5, 'Color', '#e40613')
legend('f_{H_2}','f_{H_2}^{dm}')
ylabel('H_2 [g/min]')
xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])

subplot(3,1,3);
plot(t_opt/60, 100*x_opt/M_max, 'LineWidth', 1.5, 'Color', '#32a202')
hold on
grid on
xline(t_p(1)/60, '-','t_{p1}', 'LabelVerticalAlignment', 'middle')
xline(t_p(2)/60, '-','t_{p2}', 'LabelVerticalAlignment', 'middle')
legend('m_{H_2}')
ylabel('H_2 [%]')
xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])
xlabel('Time [h]')
