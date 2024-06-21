%% Day 5
%% Setup
syms theta1(t) theta2(t) theta3(t) x(t) m0 m1 m2 m3 l1 l2 l3 g
q0_d = 0;
q1_d = 0;
q2_d = 0;
q3_d = 0;

Kp = 10;
Ki = 300;
Kd = 20;

%% Moving Base
%% Dynamics of Motion
tau = [m0*diff(x(t), t, t) - (m2*(2*l1*cos(theta1(t))*diff(theta1(t), t)^2 + l2*sin(theta1(t) + theta2(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t)) - 2*diff(x(t), t, t) + 2*l1*sin(theta1(t))*diff(theta1(t), t, t) + l2*cos(theta1(t) + theta2(t))*(diff(theta1(t), t) + diff(theta2(t), t))^2))/2 - (m1*(l1*cos(theta1(t))*diff(theta1(t), t)^2 - 2*diff(x(t), t, t) + l1*sin(theta1(t))*diff(theta1(t), t, t)))/2 - (m3*(2*l1*cos(theta1(t))*diff(theta1(t), t)^2 + 2*l2*sin(theta1(t) + theta2(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t)) - 2*diff(x(t), t, t) + l3*cos(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t) + diff(theta2(t), t) + diff(theta3(t), t))^2 + l3*sin(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t) + diff(theta3(t), t, t)) + 2*l1*sin(theta1(t))*diff(theta1(t), t, t) + 2*l2*cos(theta1(t) + theta2(t))*(diff(theta1(t), t) + diff(theta2(t), t))^2))/2;
       g*m2*(l1*cos(theta1(t)) + (l2*cos(theta1(t) + theta2(t)))/2) + m3*(l1*sin(theta1(t)) + (l3*sin(theta1(t) + theta2(t) + theta3(t)))/2 + l2*sin(theta1(t) + theta2(t)))*(l1*cos(theta1(t))*diff(theta1(t), t)^2 + l2*sin(theta1(t) + theta2(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t)) - diff(x(t), t, t) + (l3*cos(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t) + diff(theta2(t), t) + diff(theta3(t), t))^2)/2 + (l3*sin(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t) + diff(theta3(t), t, t)))/2 + l1*sin(theta1(t))*diff(theta1(t), t, t) + l2*cos(theta1(t) + theta2(t))*(diff(theta1(t), t) + diff(theta2(t), t))^2) + (l2^2*m2*(2*diff(theta1(t), t, t) + 2*diff(theta2(t), t, t)))/24 + m2*(l1*sin(theta1(t)) + (l2*sin(theta1(t) + theta2(t)))/2)*(l1*cos(theta1(t))*diff(theta1(t), t)^2 + (l2*sin(theta1(t) + theta2(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t)))/2 - diff(x(t), t, t) + l1*sin(theta1(t))*diff(theta1(t), t, t) + (l2*cos(theta1(t) + theta2(t))*(diff(theta1(t), t) + diff(theta2(t), t))^2)/2) + g*m3*(l1*cos(theta1(t)) + (l3*cos(theta1(t) + theta2(t) + theta3(t)))/2 + l2*cos(theta1(t) + theta2(t))) + m3*(l1*cos(theta1(t)) + (l3*cos(theta1(t) + theta2(t) + theta3(t)))/2 + l2*cos(theta1(t) + theta2(t)))*(l2*cos(theta1(t) + theta2(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t)) + (l3*cos(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t) + diff(theta3(t), t, t)))/2 - l1*sin(theta1(t))*diff(theta1(t), t)^2 + l1*cos(theta1(t))*diff(theta1(t), t, t) - (l3*sin(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t) + diff(theta2(t), t) + diff(theta3(t), t))^2)/2 - l2*sin(theta1(t) + theta2(t))*(diff(theta1(t), t) + diff(theta2(t), t))^2) + m2*(l1*cos(theta1(t)) + (l2*cos(theta1(t) + theta2(t)))/2)*((l2*cos(theta1(t) + theta2(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t)))/2 - l1*sin(theta1(t))*diff(theta1(t), t)^2 + l1*cos(theta1(t))*diff(theta1(t), t, t) - (l2*sin(theta1(t) + theta2(t))*(diff(theta1(t), t) + diff(theta2(t), t))^2)/2) + (l3^2*m3*(2*diff(theta1(t), t, t) + 2*diff(theta2(t), t, t) + 2*diff(theta3(t), t, t)))/24 + (l1^2*m1*diff(theta1(t), t, t))/12 + (l1*m1*sin(theta1(t))*((l1*cos(theta1(t))*diff(theta1(t), t)^2)/2 - diff(x(t), t, t) + (l1*sin(theta1(t))*diff(theta1(t), t, t))/2))/2 + (l1^2*m1*cos(theta1(t))^2*diff(theta1(t), t, t))/4 + (g*l1*m1*cos(theta1(t)))/2 - (l1^2*m1*cos(theta1(t))*sin(theta1(t))*diff(theta1(t), t)^2)/4;
      (l2^2*m2*(2*diff(theta1(t), t, t) + 2*diff(theta2(t), t, t)))/24 + g*m3*((l3*cos(theta1(t) + theta2(t) + theta3(t)))/2 + l2*cos(theta1(t) + theta2(t))) + m3*((l3*cos(theta1(t) + theta2(t) + theta3(t)))/2 + l2*cos(theta1(t) + theta2(t)))*(l2*cos(theta1(t) + theta2(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t)) + (l3*cos(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t) + diff(theta3(t), t, t)))/2 - l1*sin(theta1(t))*diff(theta1(t), t)^2 + l1*cos(theta1(t))*diff(theta1(t), t, t) - (l3*sin(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t) + diff(theta2(t), t) + diff(theta3(t), t))^2)/2 - l2*sin(theta1(t) + theta2(t))*(diff(theta1(t), t) + diff(theta2(t), t))^2) + (l3^2*m3*(2*diff(theta1(t), t, t) + 2*diff(theta2(t), t, t) + 2*diff(theta3(t), t, t)))/24 + m3*((l3*sin(theta1(t) + theta2(t) + theta3(t)))/2 + l2*sin(theta1(t) + theta2(t)))*(l1*cos(theta1(t))*diff(theta1(t), t)^2 + l2*sin(theta1(t) + theta2(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t)) - diff(x(t), t, t) + (l3*cos(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t) + diff(theta2(t), t) + diff(theta3(t), t))^2)/2 + (l3*sin(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t) + diff(theta3(t), t, t)))/2 + l1*sin(theta1(t))*diff(theta1(t), t, t) + l2*cos(theta1(t) + theta2(t))*(diff(theta1(t), t) + diff(theta2(t), t))^2) + (g*l2*m2*cos(theta1(t) + theta2(t)))/2 + (l2*m2*sin(theta1(t) + theta2(t))*(l1*cos(theta1(t))*diff(theta1(t), t)^2 + (l2*sin(theta1(t) + theta2(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t)))/2 - diff(x(t), t, t) + l1*sin(theta1(t))*diff(theta1(t), t, t) + (l2*cos(theta1(t) + theta2(t))*(diff(theta1(t), t) + diff(theta2(t), t))^2)/2))/2 + (l2*m2*cos(theta1(t) + theta2(t))*((l2*cos(theta1(t) + theta2(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t)))/2 - l1*sin(theta1(t))*diff(theta1(t), t)^2 + l1*cos(theta1(t))*diff(theta1(t), t, t) - (l2*sin(theta1(t) + theta2(t))*(diff(theta1(t), t) + diff(theta2(t), t))^2)/2))/2;
      (l3^2*m3*(2*diff(theta1(t), t, t) + 2*diff(theta2(t), t, t) + 2*diff(theta3(t), t, t)))/24 + (g*l3*m3*cos(theta1(t) + theta2(t) + theta3(t)))/2 + (l3*m3*cos(theta1(t) + theta2(t) + theta3(t))*(l2*cos(theta1(t) + theta2(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t)) + (l3*cos(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t) + diff(theta3(t), t, t)))/2 - l1*sin(theta1(t))*diff(theta1(t), t)^2 + l1*cos(theta1(t))*diff(theta1(t), t, t) - (l3*sin(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t) + diff(theta2(t), t) + diff(theta3(t), t))^2)/2 - l2*sin(theta1(t) + theta2(t))*(diff(theta1(t), t) + diff(theta2(t), t))^2))/2 + (l3*m3*sin(theta1(t) + theta2(t) + theta3(t))*(l1*cos(theta1(t))*diff(theta1(t), t)^2 + l2*sin(theta1(t) + theta2(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t)) - diff(x(t), t, t) + (l3*cos(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t) + diff(theta2(t), t) + diff(theta3(t), t))^2)/2 + (l3*sin(theta1(t) + theta2(t) + theta3(t))*(diff(theta1(t), t, t) + diff(theta2(t), t, t) + diff(theta3(t), t, t)))/2 + l1*sin(theta1(t))*diff(theta1(t), t, t) + l2*cos(theta1(t) + theta2(t))*(diff(theta1(t), t) + diff(theta2(t), t))^2))/2];

M = [diff(tau, diff(x, t, t)), ...
     diff(tau, diff(theta1, t, t)), ...
     diff(tau, diff(theta2, t, t)), ...
     diff(tau, diff(theta3, t, t))];
M = simplify(M);

gravity_term = diff(tau, g)*g;
gravity_term = simplify(gravity_term);

h = tau - M*[diff(x(t), t, t); diff(theta1(t), t, t); diff(theta2(t), t, t); diff(theta3(t), t, t)] - gravity_term;
h = simplify(h);
%% Displaying the Result
figure(1)
p = plot(t, q0, t, q1, t, q2, t, q3);
legend('Base', 'Link 1', 'Link 2', 'Link 3', 'Location', 'best')


% plot(t, q3, 'Color', [0.4940 0.1840 0.5560])
% hold on
% plot(t, q2, 'Color', [0.9290 0.6940 0.1250])
% plot(t, q1, 'Color', [0.8500 0.3250 0.0980])
% plot(t, q0, 'Color', [0 0.4470 0.7410]);
% hold off
% legend('Link 3', 'Link 2', 'Link 1', 'Base', 'Location', 'best')

xlabel('Time [s]')
ylabel('Angular Displacement [rads]')
title('Time vs. Joint Parameters')

% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(3).LineWidth = 2;
% title('PID + Computed Torque')
% set(gca,'Color','k')
