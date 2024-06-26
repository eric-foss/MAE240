%% MAE 240 Assignment 4 Problem 2
% Eric Foss
% A17068006
clear; close all; clc;


%% Problem 1

% Numerical Integration
tspan = [0 50];
x0 = [1; 1; 1; 1; 1; 1];
[t, x] = ode45(@H3BP_dynamics, tspan, x0, odeset('AbsTol',1e-13,'RelTol',1e-13));

% Jacobi Integral
J = zeros(height(x), 1);
for i = 1:height(x)

    J(i) = 0.5*(x(i, 4)^2 + x(i, 5)^2 + x(i, 6)^2) - 1/sqrt(x(i, 1)^2 + x(i, 2)^2 + x(i, 3)^2) -0.5*(3*x(i, 1)^2 - x(i, 3)^2);

end


%% Problem 2

qx = (1/3)^(1/3);
qy = 0;
qz = 0;
q = sqrt(qx^2 + qy^2 + qz^2);

Jnah = [zeros(3, 3), eye(3); -eye(3), zeros(3, 3)];
zhat = [0 -1 0; 1 0 0; 0 0 0];
zzhat = [-1 0 0; 0 -1 0; 0 0 0];

Uqq = zeros(3, 3);
Uqq(1, 1) = 3 - 1/q^3 + 3*qx^2/q^5;
Uqq(2, 2) = -1/q^3;
Uqq(3, 3) = -1 - 1/q^3;

Hxx = [-zzhat-Uqq, zhat; -zhat, eye(3)];

JHxx = Jnah*Hxx;

sigma = eig(JHxx);

[V, D] = eig(JHxx);

prop_stable = V(4:6, 1)*0.1;
prop_unstable = V(4:6, 2)*0.1;
q0 = [qx; qy; qz];

x01 = [q0; prop_stable];
x02 = [q0; -prop_stable];
x03 = [q0; prop_unstable];
x04 = [q0; -prop_unstable];

[t1, x1] = ode45(@H3BP_dynamics, flip(tspan), x01, odeset('AbsTol',1e-13,'RelTol',1e-13));
[t2, x2] = ode45(@H3BP_dynamics, flip(tspan), x02, odeset('AbsTol',1e-13,'RelTol',1e-13));
[t3, x3] = ode45(@H3BP_dynamics, tspan, x03, odeset('AbsTol',1e-13,'RelTol',1e-13));
[t4, x4] = ode45(@H3BP_dynamics, tspan, x04, odeset('AbsTol',1e-13,'RelTol',1e-13));



%% Figures
figure(1);
plot3(x(:, 1), x(:, 2), x(:, 3));
title('Hill 3BP Spacecraft Dynamics');
xlabel('Non-dimensional X'); ylabel('Non-dimensional Y'); zlabel('Non-dimensional Z');

figure(2);
plot(t, J/J(1));
title('Jacobi Integral');
ylabel('Normalized Jacobi'); xlabel('Non-dimensional Time');

figure(3); hold all;

plot3(x1(:, 1), x1(:, 2), x1(:, 3));
plot3(x2(:, 1), x2(:, 2), x2(:, 3));
plot3(x3(:, 1), x3(:, 2), x3(:, 3));
plot3(x4(:, 1), x4(:, 2), x4(:, 3));
title('Manifold Directions');
legend('Positive Stable', 'Negative Stable', 'Positive Unstable', 'Negative Unstable');
xlabel('Non-dimensional X'); ylabel('Non-dimensional Y'); zlabel('Non-dimensional Z');

