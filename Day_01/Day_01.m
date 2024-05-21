%% Day 1 
%% Pendulum with Fixed Base
figure(1)
plot(t, q1, t, q2, t, q3)
title('Time vs. Joint Parameters')
legend('Link 1', 'Link 2', 'Link 3', 'Location', 'best')
xlabel('Time [s]')
ylabel('Angular Displacement [rads]')

%% Pendulum with Moving Base
figure(1)
plot(t, x*10, t, q1, t, q2, t, q3)
title('Time vs. Joint Parameters')
legend('Base * 10', 'Link 1', 'Link 2', 'Link 3', 'Location', 'best')
xlabel('Time [s]')
ylabel('Angular Displacement [rads]')