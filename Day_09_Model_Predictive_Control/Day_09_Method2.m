%% Day 9
%% Setup
% maxInitAngle = 45 degree
% x trajectory should not become higher than lim_x
num_state           = 2;

Trajectory          = zeros(num_state, 20000);

t = (0:19999) * T_s;
amplitude = 0.5;
frequency = 0.1;
centerValue = 0;
Trajectory(1, :) = amplitude * sin(2 * pi * frequency * t) + centerValue;

L0                  = 0.01;
W0                  = 0.01;
T0                  = 0.01;
D0                  = 2700;
M0                  = L0*W0*T0*D0;

L1                  = 0.14;
W1                  = 0.01;
T1                  = 0.01;
D1                  = 2700;
M1                  = L1*W1*T1*D1;

lim_del_x           = inf;
lim_del_theta       = inf;
lim_del_dot_x       = inf;
lim_del_dot_theta   = inf;
lim_x               = 1;
lim_theta           = pi;
lim_dot_x           = inf;
lim_dot_theta       = inf;

lim_u               = 1;

disturbance         = 0;

T_s                 = 0.01;
N_p                 = 50;
M_c                 = 5;

Q           = [1 0; 0 1];
R           = 1;

Q_Bar       = kron(eye(N_p), Q);
R_Bar       = kron(eye(M_c), R);

%%
% mpcobj = mpc(plant,Ts,N_p,M_c);

%% Equation of Motion
syms x(t) theta(t) m0 m1 l1 g F
I1 = 1/12*m1*l1^2;
% Base Velocity
x_dot = diff(x);

% Joint's Angular Velocity
theta1_dot = formula(diff(theta));

S1 = sin(theta);   C1 = cos(theta);

% Linear Displacement
s0 = [x, 0].';
s1 = [x - 1/2*l1*S1, 1/2*l1*C1].';

% Linear Velocity
s0_dot = diff(s0);
s1_dot = diff(s1);

% Kinematic Energy
KE0 = 1/2*m0*(s0_dot.')*s0_dot;
KE1 = 1/2*m1*(s1_dot.')*s1_dot + 1/2*I1*theta1_dot^2;

% Potential Energy
PE0 = 0;
PE1 = 1/2*m1*g*l1*C1;

% Lagrangian
L = KE0 + KE1 - PE0 - PE1;

dL_dx_dot = diff(L, x_dot);
diff_dL_dx_dot = diff(dL_dx_dot, t);
dL_dx = diff(L, x);

dL_dtheta1_dot = diff(L, theta1_dot);
diff_dL_dtheta1_dot = diff(dL_dtheta1_dot, t);
dL_dtheta1 = diff(L, theta);

tau = [diff_dL_dx_dot - dL_dx;
       diff_dL_dtheta1_dot - dL_dtheta1];
tau = simplify(tau);

M = [diff(tau, diff(x, t, t)), ...
     diff(tau, diff(theta, t, t))];
M = simplify(M);

gravity_term = diff(tau, g)*g;
gravity_term = simplify(gravity_term);

h = tau - M*[diff(x(t), t, t); diff(theta(t), t, t)] - gravity_term;
h = simplify(h);

acc = simplify(M^-1*([F; 0] - gravity_term - h));

g_x = simplify(diff(acc, F)*F);
f_x = simplify(acc - g_x);

f_x  = [diff(x, t);
        diff(theta, t);
        f_x];

g_x = [0;
       0;
       g_x];

%% Euler Discretization
A_d = formula((eye(4) + T_s*A_E));
B_d = formula(T_s*B_E);
C_d = sym([eye(2), zeros(2)]);

%% Augmented State Space
A_a = formula([A_d,     zeros(4, 2);
               C_d*A_d, eye(2, 2)]);
B_a = formula([B_d; C_d*B_d]);
C_a = sym([zeros(2, 4), eye(2)]);

%% Prediction Model
[m_A, n_A] = size(A_a);
[m_C, n_C] = size(C_a);

F_a = sym(zeros(N_p*m_C, n_A));
h_a = sym(zeros(N_p*m_C, n_A));

F_a(1:2, 1:6) = formula(C_a*A_a);
h_a(1:2, 1:6) = formula(C_a);

for i = 2 : N_p
    F_a(i*2-1:i*2, :) = formula(F_a((i-1)*2-1:(i-1)*2,:)*A_a);
    h_a(i*2-1:i*2, :) = formula(h_a((i-1)*2-1:(i-1)*2,:)*A_a);
end

v = h_a* B_a;
Phi = v;

for i = 2 : M_c
    Phi(:, i) = [zeros(2*(i-1), 1); v(1:2*(N_p-i+1), 1)];
end

%% Simulation
x_0 = [0, pi/20, 0, 0];
ref = reshape(zeros(2, N_p), [], 1);

for t = T_s:T_s:10
    
end




%% Displaying the Result
% Output y
figure(1)
clf;
yyaxis left
h2 = plot(t, x2, '-b', 'DisplayName', 'theta');
hold on
ylabel('\theta [rad]')
set(gca, 'YColor', 'b');
yyaxis right
h1 = plot(t, x1, '-r', 'DisplayName', 'x');
ylabel('x [m]')
set(gca, 'YColor', 'r');
xlabel('t [s]')
title('State vs. Time')
legend([h1, h2], {'x', '\theta'}, 'Location', 'best');
hold off

% Sliding Surface
figure(2)
fplot(@(t) -C_x*t)
hold on
quiver(x1(1:end-1), x3(1:end-1), diff(x1), diff(x3), 1, 'r', 'LineWidth', 1);
plot(0, 0, 'ok', 'MarkerSize', 20)
xlim([min(x1)-0.01, max(x1)]*1.1)
ylim([min(x3), max(x3)]*1.1)
title('Phase Portrait: Linear Motion')
xlabel('x1 [m]')
ylabel('x3 [m/s]')
legend('s = 0', 'phase portrait', '', 'Location', 'best')
hold off

figure(3)
plot(C_theta*(x2-z), x4)
hold on
quiver(x2(1:end-1), x4(1:end-1), diff(x2), diff(x4), 1, 'r', 'LineWidth', 1);
plot(0, 0, 'ok', 'MarkerSize', 20)
% xlim([-1, 1])
title('Phase Portrait: Angular Motion')
xlabel('x2 [rad]')
ylabel('x4 [rad/s]')
legend('s = 0', 'phase portrait', '', 'Location', 'best')
hold off

figure(4)
plot(t, S_theta)
hold on
plot(t, S_x)
title('The magnitude of Sliding Surface')
xlabel('t [s]')
ylabel('S')
legend('S_{\theta}', 'S_x', 'Location', 'best')
hold off

%% Comparing Parameters
figure(1)
plot(t, x1);
hold on

figure(2)
plot(t, x2);
hold on

%% Legends
figure(1)
title('Effect of C_x on x')
xlabel('t [s]')
ylabel('x [m]')
hold off
legend('C_x = 0.1', 'C_x = 1', 'C_x = 5', 'C_x = 10', 'Location', 'best')

figure(2)
title('Effect of C_x on \theta')
xlabel('t [s]')
ylabel('\theta [rads]')
hold off
legend('C_x = 0.1', 'C_x = 1', 'C_x = 5', 'C_x = 10', 'Location', 'best')


%% Input and Output
figure(1)
clf;
yyaxis left
h2 = plot(t, x2, '-b');
hold on
h1 = plot(t, x1, '-r', 'DisplayName', 'x');
ylabel('Ouput x and \theta')
set(gca, 'YColor', 'b');
yyaxis right
h3 = plot(t, u, '-g');
ylabel('Input F')
set(gca, 'YColor', 'g');
xlabel('t [s]')
title('State vs. Time')
legend([h1, h2, h3], {'x', '\theta', 'u'}, 'Location', 'best');
hold off

%% Reference and Output
figure(1)
plot(t, ref, '-b');
hold on
plot(t, x1, '-r');
ylabel('x [m]')
xlabel('t [s]')
title('Sinusoidal Response')
legend('reference', 'x', 'Location', 'best');
hold off


