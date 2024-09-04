%% Day 9 - Method 1
%% Setup
% maxInitAngle = 45 degree
% x trajectory should not become higher than lim_x
T_s                 = 0.01;
num_state           = 2;

Trajectory          = zeros(num_state, 20000);
% Trajectory          = load("150deg.mat").temp';
t = (0:19999) * T_s;
amplitude           = 0.3;
frequency           = 0.4;
centerValue         = 0;
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

lim_x               = 1;
lim_theta           = 2*pi;

lim_u               = 1;

disturbance         = 0; % [3.43, 3.35, 2.92, 2.8, 0.79, 1.3, 0.98]

N_p                 = 50; % [23, 30, 40, 50, 80, 60, 70]
M_c                 = 5;

Q                   = [1 0; 0 1];
R                   = 1; 

Q_Bar               = kron(eye(N_p), Q);
R_Bar               = kron(eye(M_c), R);

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
A_E = formula(jacobian(f_x, [x, theta, diff(x, t), diff(theta, t)]));
B_E = formula(jacobian(g_x, F));

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

%% Quadratic Programming
% Try this with very small N_p and N_c
syms del_x del_theta del_dot_x del_dot_theta
x_a = [del_x; del_theta; del_dot_x; del_dot_theta; x; theta];
R_ref = sym(zeros(N_p*2, 1));
Delta_U = formula(-(Phi.'*Q_Bar*Phi + R_Bar)^(-1)*Phi.'*Q_Bar*(F_a*x_a - R_ref));
delta_u_k = Delta_U(1);

%% Comparing Parameters
figure(1)
plot(t, x1);
hold on

figure(2)
plot(t, x2);
hold on

figure(3)
plot(t(1:end-1), computeCumulativeWork(x1, u));
hold on

%% Legends
figure(1)
title('Effect of R on x')
xlabel('t [s]')
ylabel('x [m]')
hold off
legend('M = 2', 'M = 5', 'M = 15', 'M = 45', 'Location', 'best')

figure(2)
title('Effect of R on \theta')
xlabel('t [s]')
ylabel('\theta [rads]')
hold off
legend('M = 2', 'M = 5', 'M = 15', 'M = 45', 'Location', 'best')


figure(3)
title('Effect of R on Energy Input')
xlabel('t [s]')
ylabel('u [N]')
hold off
legend('M = 2', 'M = 5', 'M = 15', 'M = 45', 'Location', 'best')

%% Robustness
figure(4)
variable = [23, 30, 40, 50, 60, 70, 80];
Dist =  [3.43, 3.35, 2.92, 2.8, 1.3, 0.98, 0.79];
plot(variable, Dist, '-ob');
xlabel('R')
ylabel('Disturbance [N]')
title('Prediction Horizon (N)')

%%
computeCumulativeWork(x1, u);
ans(end)
%% Work Input
figure(4)
variable = 2:20;
Dist =  [0.0462, 0.0431, 0.0352, 0.0342, 0.0361, 0.0388, 0.0412, 0.0436, 0.0455, ...
         0.0468, 0.0468, 0.0461, 0.0452, 0.0446, 0.0442, 0.0441, 0.0441, 0.0441, 0.0441];
plot(variable, Dist, '-ob');
xlabel('Control Horizon (M)')
ylabel('Total Work [J]')

%% System Identification (Using FFT)
% Perform FFT on both signals
n = length(x1(1000:end));
fft_input = fft(Trajectory(1, 1000:length(x1)));
fft_output = fft(x1(1000:end));

% Frequency axis
fs = 1/T_s; % Sampling frequency
f = (0:n-1)*(fs/n); % Frequency vector

% Identify the frequency of interest (e.g., fundamental frequency)
[~, idx] = min(abs(f - frequency)); % Find index of the target frequency

% Calculate gain and phase shift at the target frequency
gain = abs(fft_output(idx)) / abs(fft_input(idx));
phase_shift = angle(fft_output(idx)) - angle(fft_input(idx));

% Convert phase shift to degrees
phase_shift_degrees = rad2deg(phase_shift);

% Display results
fprintf('Gain: %.4f\n', gain);
fprintf('Phase Shift: %.4f degrees\n', phase_shift_degrees);

figure(5)
plot(t(1:length(x1)), Trajectory(1, 1:length(x1)), t(1:length(x1)), x1)
xlabel('Time [s]')
ylabel('Position [m]')
legend('x_d', 'x', 'location', 'best')

%% Bode Plot
freqs = [0.01, 0.05, 0.1, ...
         0.2, 0.3, 0.4, ...
         0.5, 0.8, 1, 5];
gains = [1, 1.02, 1.06, ...
         1.2736, 1.5181, 1.1383, ...
         0.7023, 0.25889, 0.1706, 0.0036];
phases = [-0.03, -0.03, -1.14, ...
          -8.4141, -33.5539, -66.7598, ... 
          -78.2385, -66.7061, -50.1954, -205.5448];

% Convert gain to dB
gains_dB = 20 * log10(gains);

figure(6)
% Plot Gain (Magnitude) in dB
subplot(2,1,1);
semilogx(freqs, gains_dB, 'b-o');
grid on;
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
xlim([min(freqs) max(freqs)]); 
% Plot Phase in degrees
subplot(2,1,2);
semilogx(freqs, phases, 'r-o');
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
xlim([min(freqs) max(freqs)]); 
% Improve layout
sgtitle('Bode Plot');
%% Helper Function
function cumulativeWork = computeCumulativeWork(x, F)
    % Function to compute the cumulative sum of work done by the input force
    %
    % Inputs:
    %   F - Vector of input forces
    %   x - Vector of cart positions corresponding to the forces
    %
    % Output:
    %   cumulativeWork - The cumulative sum of work done by the force F

    % Calculate the change in position (displacement)
    delta_x = diff(x);

    % Calculate the work done in each segment
    work = abs(F(1:end-1) .* delta_x);

    % Calculate the cumulative sum of work done
    cumulativeWork = cumsum(work);
end
