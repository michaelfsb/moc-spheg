%% Single-phase

import casadi.*

clearvars
load("inputs\data.mat")
addpath('models')

%% Problem constants

Tf = 1440;   % Final time (min)
N = 42;      % Number of control intervals

M_0 = 164;     % Initial volume of hydrogen (g)
M_min = 49.2;  % Minimum volume of hydrogen (g)
M_max = 205;   % Maximum volume of hydrogen (g)
I_e_0 = 30;    % Initial current (A)
I_e_min = 30;  % Minimum current (A)
I_e_max = 100; % Maximum current (A)

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
v_h2_dot = f_h2 - HydrogenDemand(time);

% Lagrange cost function
f_l = (p_el - p_ps)^2;

% Right-hand side
f = casadi.Function('f', ...
    {v_h2, i_el, time}, {v_h2_dot, f_l}, ...
    {'x', 'u', 't'}, {'x_dot', 'L'});

%% Create the NPL

% Time grid
h = Tf/N;
nGrid = 2*N-1;
t = linspace(0,Tf,N);

X = casadi.MX.sym('X',nGrid);
U = casadi.MX.sym('U',nGrid);

L = 0;
g = MX(nGrid-1,1);
for k=1:2:nGrid-1
    i = ceil(k/2);
    [f_k_0, w_k_0] = f(X(k), U(k), t(i));
    [f_k_1, w_k_1] = f(X(k+1), U(k+1), t(i)+h/2);
    [f_k_2, w_k_2] = f(X(k+2), U(k+2), t(i+1));
    
    L = L + h*(w_k_0 + 4*w_k_1 + w_k_2)/6;
    
    g(k) =  X(k+2) - X(k) - h*(f_k_0 + 4*f_k_1 + f_k_2)/6;
    g(k+1) =  X(k+1) - (X(k+2) + X(k))/2 - h*(f_k_0 - f_k_2)/8;
end
lbg = zeros(length(g),1);
ubg = zeros(length(g),1);

%% Set Path constrains and initial conditions

x = [X;U];

lbx = [M_min*ones(nGrid,1); ...
    I_e_min*ones(nGrid,1)];
ubx = [M_max*ones(nGrid,1); ...
    I_e_max*ones(nGrid,1)];

% Set the initial condition for the state
lbx(1) = M_0;
ubx(1) = M_0;

% Set guess
xGuess = (M_max - M_min)/2;
uGuess = (I_e_max - I_e_min)/2;
x0 = [xGuess*ones(length(X),1) ; uGuess*ones(length(U),1)];

%% Solve

nlp = struct('x', x, 'f', L, 'g', g);
solver = nlpsol('solver', 'ipopt', nlp);

solution = solver('x0', x0, 'lbx', lbx, 'ubx', ubx,...
    'lbg', lbg, 'ubg', ubg);

%% Get solution

w_opt = (full(solution.x));
x_opt = w_opt(1:end/2);
u_opt = w_opt(end/2+1:end);
t_opt = linspace(0,1440, length(u_opt));

%% Simulate solution

nOpt = length(u_opt);
f_h2_prd_opt =  zeros(nOpt,1);
v_el_opt =  zeros(nOpt,1);
p_el_opt =  zeros(nOpt,1);
p_ps_opt = zeros(nOpt,1);
i_ps_opt = zeros(nOpt,1);
for i=1:nOpt
    [f_h2_prd_opt(i), v_el_opt(i), p_el_opt(i)] = electrolyzer_model(u_opt(i));
    [i_ps_opt(i), ~, p_ps_opt(i)] = photovoltaic_model(full(Irradiation(t_opt(i))));
end

fprintf('\nCost: %f W^2\n',full(solution.f))
fprintf('Grid energy: %4.2f W\n',sqrt(full(solution.f)))

%% Plots

figure()
subplot(3,1,1);
plot(t_opt/60, p_el_opt, 'LineWidth', 1.5, 'Color', '#0152a1')
hold on
grid on
plot(t_opt/60, p_ps_opt, 'LineWidth', 1.5, 'Color', '#ff8700')
legend('P_E','P_S')
ylabel('Power [W]')
xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])

subplot(3,1,2);
plot(t_opt/60, f_h2_prd_opt, 'LineWidth', 1.5, 'Color', '#0152a1')
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
legend('m_{H_2}')
ylabel('H_2 [%]')
xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])
xlabel('Time [h]')