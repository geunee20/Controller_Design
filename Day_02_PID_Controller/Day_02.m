%% Day 2
%% PID Controller with Fixed Base
%% Setup
qb_d = 0;
q1_d = 0;
q2_d = 0;
q3_d = 0;

%% Displaying the Result
figure(1)
plot(t, q1, t, q2, t, q3)
title('Time vs. Joint Parameters')
legend('Link 1', 'Link 2', 'Link 3', 'Location', 'best')
xlabel('Time [s]')
ylabel('Angular Displacement [rads]')

%% PID Controller with Moving Base
%% Displaying the Result
figure(1)
plot(t, qb, t, q1, t, q2, t, q3)
title('Time vs. Joint Parameters')
legend('Base', 'Link 1', 'Link 2', 'Link 3', 'Location', 'best')
xlabel('Time [s]')
ylabel('Angular Displacement [rads]')