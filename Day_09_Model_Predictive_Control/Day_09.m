%% Day 9
%% Setup
num_state   = 4;

Trajectory  = zeros(num_state, 20001);

L0          = 0.01;
W0          = 0.01;
T0          = 0.01;
D0          = 2700;
M0          = L0*W0*T0*D0;

L1          = 0.14;
W1          = 0.01;
T1          = 0.01;
D1          = 2700;
M1          = L1*W1*T1*D1;

limit       = 0.45;

disturbance = 0;

T_s         = 0.001;
N_p         = 30;
M_c         = 15;

Q           = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
R           = 0.01;

Q_Bar       = createBlockDiagonal(Q, N_p);
R_Bar       = createBlockDiagonal(R, M_c);
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

%% Jacobian Linearization
A_E = jacobian(f_x, [x, theta, diff(x, t), diff(theta, t)]);
B_E = jacobian(g_x, F);

%% Augmented State Space


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


function bar_Mat = createBlockDiagonal(Q, N)
    % createBlockDiagonal - Creates a block diagonal matrix with Q repeated N times
    %
    % Syntax: bar_Q = createBlockDiagonal(Q, N)
    %
    % Inputs:
    %    Q - The matrix to be repeated on the diagonal
    %    N - The number of times to repeat Q on the diagonal
    %
    % Outputs:
    %    bar_Mat - The resulting block diagonal matrix with Q repeated N times

    % Create a cell array with N copies of Q
    Q_cell = repmat({Q}, 1, N);

    % Create the block diagonal matrix \bar{Q}
    bar_Mat = blkdiag(Q_cell{:});
end
