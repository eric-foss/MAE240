%% MAE 240 Assignment 3 Problem 1
% Eric Foss
% A17068006

const.nu = 0.0122;


r0 = [10000; 0; 0];
v0 = [0; 400; 0];

x0 = [r0'; v0'];

spy = 3.154e+7;

time = [0, 0.1*spy];


[t, x] = ode45(@CR3BP_dynamics, time , x0, odeset('AbsTol',1e-11,'RelTol',1e-9), const);

%plot3(x(1), x(2), x(3))