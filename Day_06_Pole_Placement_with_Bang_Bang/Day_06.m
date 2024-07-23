%% Day 6
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

Disturbance = 0;

limit = 0.45;

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
dL_dtheta1 = diff(L, theta1);

tau = [diff_dL_dx_dot - dL_dx;
       diff_dL_dtheta1_dot - dL_dtheta1];
tau = simplify(tau);
tau_linear = subs(tau, [cos(theta1), sin(theta1), diff(theta1, t)^2], [1, theta1, 0]);

%% Linearization
M = [diff(tau_linear, diff(x, t, t)), ...
     diff(tau_linear, diff(theta1, t, t))];
M = simplify(M);

gravity_term = diff(tau_linear, g)*g;
gravity_term = simplify(gravity_term);

acc = simplify(M^-1*([F; 0] - gravity_term));

acc = subs(acc, [m0, m1, l1], [M0, M1, L1]);
A = diff(acc, theta1(t));
A = double(subs(A, g, 9.81));
A = [0 0 1 0;
     0 0 0 1;
     0 A(1) 0 0;
     0 A(2) 0 0];
B = double(diff(acc, F));
B = [0; 0; B(1); B(2)];
C = [1 0 0 0;
     0 1 0 0];
D = 0;

%% Controllability
co = ctrb(A, B);
controllability = rank(co);

%% Gain matrix from Pole Placement Method
K_place = place(A, B, [-15.000, ... 
                       -15.001, ...
                       -15.002, ...
                       -15.003]);

%% Displaying the Result 2
figure(1)
h1 = plot(t, x1);
hold on

figure(2)
h2 = plot(t, x2);
hold on

%%
figure(1)
legend('s=-1', 's=-5', 's=-10', 's=-15', 'location', 'best')
title('x_1 vs. time')
xlabel('time [s]')
ylabel('x_1 [m]')
figure(2)
title('Linear Displacement')
legend('s=-5', 's=-10', 's=-15', 's=-20', 'location', 'best')
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
