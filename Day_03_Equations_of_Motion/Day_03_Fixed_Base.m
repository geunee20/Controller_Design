%% variablees
syms theta1(t) theta2(t) theta3(t) l1 l2 l3 m1 m2 m3 g

I1 = 1/12*m1*l1^2;
I2 = 1/12*m2*l2^2;
I3 = 1/12*m3*l3^2;

%% Kinematics
% Joint's Angular Velocity
theta1_dot = diff(theta1);
theta2_dot = diff(theta2);
theta3_dot = diff(theta3);

% Body's Angular Velocity
omega1 = theta1_dot;
omega2 = theta1_dot + theta2_dot;
omega3 = theta1_dot + theta2_dot + theta3_dot;

S1 = sin(theta1);   S12 = sin(theta1 + theta2);   S123 = sin(theta1 + theta2 + theta3);
C1 = cos(theta1);   C12 = cos(theta1 + theta2);   C123 = cos(theta1 + theta2 + theta3);

% Linear Displacement
s1 = [1/2*l1*C1, 1/2*l1*S1].';
s2 = [l1*C1 + 1/2*l2*C12, l1*S1 + 1/2*l2*S12].';
s3 = [l1*C1 + l2*C12 + 1/2*l3*C123, l1*S1 + l2*S12 + 1/2*l3*S123].';

% Linear Velocity
s1_dot = diff(s1);
s2_dot = diff(s2);
s3_dot = diff(s3);

%% Energy
% Kinematic Energy
KE1 = 1/2*m1*(s1_dot.')*s1_dot + 1/2*I1*omega1^2;
KE2 = 1/2*m2*(s2_dot.')*s2_dot + 1/2*I2*omega2^2;
KE3 = 1/2*m3*(s3_dot.')*s3_dot + 1/2*I3*omega3^2;

% Potential Energy

PE1 = m1*g*1/2*l1*S1;
PE2 = m2*g*(l1*S1 + 1/2*l2*S12);
PE3 = m3*g*(l1*S1 + l2*S12 + 1/2*l3*S123);

L = KE1 + KE2 + KE3 - PE1 - PE2 - PE3;

dL_dtheta1_dot = diff(L, theta1_dot);
diff_dL_dtheta1_dot = diff(dL_dtheta1_dot, t);
dL_dtheta1 = diff(L, theta1);

dL_dtheta2_dot = diff(L, theta2_dot);
diff_dL_dtheta2_dot = diff(dL_dtheta2_dot, t);
dL_dtheta2 = diff(L, theta2);

dL_dtheta3_dot = diff(L, theta3_dot);
diff_dL_dtheta3_dot = diff(dL_dtheta3_dot, t);
dL_dtheta3 = diff(L, theta3);

tau = [diff_dL_dtheta1_dot - dL_dtheta1;
       diff_dL_dtheta2_dot - dL_dtheta2;
       diff_dL_dtheta3_dot - dL_dtheta3];

%% Jacobian
J = jacobian([s3; theta1+theta2+theta3], [theta1, theta2, theta3]);