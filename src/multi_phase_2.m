%% Multi-phase
import casadi.*

clearvars
addpath('inputs')
addpath('models')
addpath('interpolation')
load("data.mat")

%% Problem constants

Tf = 1440;   % Final time (min)
N = 14;      % Number of control intervals per phase

M_0 = 164;     % Initial volume of hydrogen (g)
M_min = 49.2;  % Minimum volume of hydrogen (g)
M_max = 205;   % Maximum volume of hydrogen (g)
I_e_0 = 30;    % Initial current (A)
I_e_min = 20;  % Minimum current (A)
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
tic
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
toc
%% Solve

nlp = struct('x', x, 'f', L, 'g', g);
solver = nlpsol('solver', 'ipopt', nlp);
solution = solver('x0', x0, 'lbx', lbx, 'ubx', ubx,...
    'lbg', lbg, 'ubg', ubg);

%% Get solution

w_sol = (full(solution.x));
t_p = w_sol(end-1:end);

t_sol = [ linspace(0,t_p(1), nGrid); ...
    linspace(t_p(1),t_p(2), nGrid); ...
    linspace(t_p(2),1440, nGrid) ];

x_sol = [ w_sol(1:nGrid)'; ...
    w_sol(nGrid+1:2*nGrid)'; ...
    w_sol(2*nGrid+1:3*nGrid)' ];

u_sol = [ w_sol(3*nGrid+1:4*nGrid)'; ...
    w_sol(4*nGrid+1:5*nGrid)'; ...
    w_sol(5*nGrid+1:6*nGrid)' ];

%% Simulate solution

f_sol = [ full(f_st(x_sol(1,:), u_sol(1,:), t_sol(1,:))); ...
    full(f_on(x_sol(2,:), u_sol(2,:), t_sol(2,:))); ...
    full(f_st(x_sol(3,:), u_sol(3,:), t_sol(3,:))) ];

nInter = 10*nGrid;
t_inter = [ linspace(0, t_p(1), nInter); ...
    linspace(t_p(1), t_p(2), nInter); ...
    linspace(t_p(2), 1440, nInter) ];

u_inter = [ interp_ctr(t_sol(1,:), u_sol(1,:), t_inter(1,:)); ...
    interp_ctr(t_sol(2,:), u_sol(2,:), t_inter(2,:)); ...
    interp_ctr(t_sol(3,:), u_sol(3,:), t_inter(3,:)) ];

x_inter = [ interp_std(t_sol(1,:), x_sol(1,:), f_sol(1,:), t_inter(1,:)); ...
    interp_std(t_sol(2,:), x_sol(2,:), f_sol(2,:), t_inter(2,:)); ...
    interp_std(t_sol(3,:), x_sol(3,:), f_sol(3,:), t_inter(3,:)) ];

f_h2_prd_opt =  zeros(3, nInter);
v_el_opt = zeros(3, nInter);
p_el_opt = zeros(3, nInter);
p_ps_opt = zeros(3, nInter);
p_grid = zeros(3, nInter);
p_grid2 = zeros(3,nInter);
for p=1:3
    for i=1:nInter
        [f_h2_prd_opt(p, i), v_el_opt(p, i), p_el_opt(p, i)] = electrolyzer_model(u_inter(p, i));
        [~, ~, p_ps_opt(p, i)] = photovoltaic_model(full(Irradiation(t_inter(p, i))));
        if p == 1 || p == 3
            f_h2_prd_opt(p, i) = 0;
        end
        if p_el_opt(p, i) > p_ps_opt(p, i) - 1
            p_grid(p, i) = p_el_opt(p, i) -  p_ps_opt(p, i);
        else
            p_grid2(p, i) = p_ps_opt(p, i) - p_el_opt(p, i);
        end
    end
end

e_grid = trapz(t_inter(1,:)/60, p_grid(1,:))/1000 + ...
    trapz(t_inter(2,:)/60, p_grid(2,:))/1000 + ...
    trapz(t_inter(3,:)/60, p_grid(3,:))/1000;

e_grid2 = trapz(t_inter(1,:)/60, p_grid2(1,:))/1000 + ...
    trapz(t_inter(2,:)/60, p_grid2(2,:))/1000 + ...
    trapz(t_inter(3,:)/60, p_grid2(3,:))/1000;

e_el = trapz(t_inter(1,:)/60, p_el_opt(1,:))/1000 + ...
    trapz(t_inter(2,:)/60, p_el_opt(2,:))/1000 + ...
    trapz(t_inter(3,:)/60, p_el_opt(3,:))/1000;

h2_prd = trapz(t_inter(1,:), f_h2_prd_opt(1,:)) + ...
    trapz(t_inter(2,:), f_h2_prd_opt(2,:)) + ...
    trapz(t_inter(3,:), f_h2_prd_opt(3,:));

fprintf('\n\nCosumed grid energy: %4.2f kWh\n', e_grid)
fprintf('Absorbed grid energy: %4.2f kWh\n', e_grid2)
fprintf('Electrolyzer energy: %4.2f kWh\n', e_el)
fprintf('Prod H2: %4.2f g\n', h2_prd)
fprintf('Final H2: %4.2f %\n', 100*x_inter(end)/M_max)

%% Plots

figure()
subplot(4,1,1);
hold on
grid on
plot(t_inter(1,:)/60, p_ps_opt(1,:)/1000, 'LineWidth', 1.5, 'Color', '#ff8700')
plot(t_inter(2,:)/60, p_ps_opt(2,:)/1000, 'LineWidth', 1.5, 'Color', '#ff8700')
plot(t_inter(3,:)/60, p_ps_opt(3,:)/1000, 'LineWidth', 1.5, 'Color', '#ff8700')
plot(t_inter(1,:)/60, p_el_opt(1,:)/1000, 'LineWidth', 1.5, 'Color', '#0152a1')
plot(t_inter(2,:)/60, p_el_opt(2,:)/1000, 'LineWidth', 1.5, 'Color', '#0152a1')
plot(t_inter(3,:)/60, p_el_opt(3,:)/1000, 'LineWidth', 1.5, 'Color', '#0152a1')
legend('P_S', '','', 'P_E')
ylabel('Power [kW]')
xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])

subplot(4,1,2);
hold on
grid on
plot((0:1:1440)/60,full(HydrogenDemand(0:1:1440)), 'LineWidth', 1.5, 'Color', '#e40613')
plot(t_inter(1,:)/60, f_h2_prd_opt(1,:), 'LineWidth', 1.5, 'Color', '#0152a1')
plot(t_inter(2,:)/60, f_h2_prd_opt(2,:), 'LineWidth', 1.5, 'Color', '#0152a1')
plot(t_inter(3,:)/60, f_h2_prd_opt(3,:), 'LineWidth', 1.5, 'Color', '#0152a1')
legend('f_{H_2}^{dm}','f_{H_2}')
ylabel('H_2 [g/min]')
xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])

subplot(4,1,3);
hold on
grid on
plot(t_inter(1,:)/60, 100*x_inter(1,:)/M_max, 'LineWidth', 1.5, 'Color', '#32a202')
plot(t_inter(2,:)/60, 100*x_inter(2,:)/M_max, 'LineWidth', 1.5, 'Color', '#32a202')
plot(t_inter(3,:)/60, 100*x_inter(3,:)/M_max, 'LineWidth', 1.5, 'Color', '#32a202')
xline(t_p(1)/60, '-','t_{p1}', 'LabelVerticalAlignment', 'middle')
xline(t_p(2)/60, '-','t_{p2}', 'LabelVerticalAlignment', 'middle')
legend('m_{H_2}')
ylabel('H_2 [%]')
xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])

subplot(4,1,4);
hold on
grid on
plot(t_inter(1,:)/60, u_inter(1,:), 'LineWidth', 1.5, 'Color', '#0152a1')
plot(t_inter(2,:)/60, u_inter(2,:), 'LineWidth', 1.5, 'Color', '#0152a1')
plot(t_inter(3,:)/60, u_inter(3,:), 'LineWidth', 1.5, 'Color', '#0152a1')
y1 = yline(I_e_min, '--','i_{E}^{min}','LineWidth',2);
y1.LabelHorizontalAlignment = 'center';
legend('i_{E}')
ylabel('Current [A]')
xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])
xlabel('Time [h]')