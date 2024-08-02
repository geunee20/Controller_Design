%% Day 8
%% Setup
x1_d = 0;
x2_d = 0;
x3_d = 0;
x4_d = 0;

L0 = 0.01;
W0 = 0.01;
T0 = 0.01;
D0 = 2700;
M0 = L0*W0*T0*D0;

L1 = 0.14;
W1 = 0.01;
T1 = 0.01;
D1 = 2700;
M1 = L1*W1*T1*D1;

K = 100;
C_x = 1;
C_theta = 3;

lambda = 1;
epsilon = 0.5;
m = 3;

limit = .45;
Disturbance = 0;

%% Equation of Motion
syms x(t) theta1(t) m0 m1 l1 g F
I1 = 1/12*m1*l1^2;
% Base Velocity
x_dot = diff(x);

% Joint's Angular Velocity
theta1_dot = formula(diff(theta1));

S1 = sin(theta1);   C1 = cos(theta1);

% Linear Displacement
s0 = [x, 0].';
s1 = [x + 1/2*l1*S1, 1/2*l1*C1].';

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
dL_dtheta1 = diff(L, theta1);

tau = [diff_dL_dx_dot - dL_dx;
       diff_dL_dtheta1_dot - dL_dtheta1];
tau = simplify(tau);

M = [diff(tau, diff(x, t, t)), ...
     diff(tau, diff(theta1, t, t))];
M = simplify(M);

gravity_term = diff(tau, g)*g;
gravity_term = simplify(gravity_term);

h = tau - M*[diff(x(t), t, t); diff(theta1(t), t, t)] - gravity_term;
h = simplify(h);

acc = simplify(M^-1*([F; 0] - gravity_term - h));

g_x = diff(acc, F)*F;
f_x = simplify(acc - g_x);
%% Displaying the Result 2
figure(1)
h1 = plot(t, x1);
hold on

figure(2)
h2 = plot(t, x2);
hold on

%%
figure(1)
legend('F = 0.6N', 'F = 1N', 'F = 3N', 'F = 10N', 'location', 'best')
title('x_1 vs. time')
xlabel('time [s]')
ylabel('x_1 [m]')
figure(2)
legend('F = 0.6N', 'F = 1N', 'F = 3N', 'F = 10N', 'location', 'best')
title('x_2 vs. time')
xlabel('time [s]')
ylabel('x_2 [rad]')

%% Displaying the Result
figure(1)

h2 = plot(t, x2, 'b', 'DisplayName', 'x2');
hold on
h1 = plot(t, x1, 'r', 'DisplayName', 'x1');
hold off

legend([h1, h2], {'x_1', 'x_2'}, 'location', 'best');

%% Energy Calculation
function E = calculate_pendulum_energy(x_dot, theta, theta_dot, m0, m1, l, g)
    % Kinetic Energy of the cart
    T_c = 0.5 * m0 * x_dot.^2;

    % Kinetic Energy of the pendulum
    T_p = (m1*x_dot.^2)/2 + (l^2*m1*theta_dot.^2)/6 - (l*m1*cos(theta).*theta_dot.*x_dot)/2;

    % Potential Energy of the pendulum
    V_p = m1 * g * 0.5 * l * cos(theta);

    % Total Energy
    E = T_c + T_p + V_p;
end