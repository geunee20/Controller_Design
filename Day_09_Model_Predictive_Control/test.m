A = eye(4);
B = eye(4, 1);
C = eye(2, 4);

x_d = zeros(6, 1);
x = zeros(6, 1);

Q = eye(20);
R = eye(3);

x_d = zeros()

[Aa, Ba, Ca] = fcn(A, B);
[F, Phi] = fcn2(Aa, Ba, Ca, 10, 3);
u = fcn3(x_d, x, F, Phi, Q, R);

function [A_a, B_a, C_a] = fcn(A_d, B_d)

C_d = [eye(2), zeros(2)];

[m, n] = size(C_d);
[~, o] = size(B_d);

A_a = eye(m+n, m+n);
A_a(1:n ,1:n) = A_d;
A_a(n+1:n+m, 1:n) = C_d*A_d;

B_a = zeros(n+m, o);
B_a(1:n, :) = B_d;
B_a(n+1:n+m, :) = C_d*B_d;

C_a = zeros(m, n+m);
C_a(:, n+1:n+m) = eye(m, m);
end

function [F, Phi] = fcn2(A_a, B_a, C_a, N, M)

[m_A, n_A] = size(A_a);
[m_C, n_C] = size(C_a);

F = zeros(N*m_C, n_A);
h = zeros(N*m_C, n_A);

F(1:2, :) = C_a*A_a;
h(1:2, :) = C_a;

for i = 2 : N
    F(i*2-1:i*2, :) = F((i-1)*2-1:(i-1)*2,:)*A_a;
    h(i*2-1:i*2, :) = h((i-1)*2-1:(i-1)*2,:)*A_a;
end

v = h * B_a;
Phi = zeros(N*2, M);
Phi(:, 1) = v;

for i = 2 : M
    Phi(:, i) = [zeros(2*(i-1), 1); v(1:2*(N-i+1), 1)];
end
end

function u = fcn3(x_d, x, F, Phi, Q, R)

U = -(Phi'*Q*Phi + R)^(-1)*Phi'*Q*(F*x-x_d);
u = U(1);
end