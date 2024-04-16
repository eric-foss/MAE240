%% MAE 240 Assignment 1 Problem 1
% Eric Foss
% A17068006

%% Universal Constants
const.G = 6.6742e-11; %Gravitation Constant [N m^2/kg^2]

%% Part A

n = input("How many bodies?");

const.m_all = rand(1, n) * 10^29; %random mass for each body [kg]

r0_all = (rand(1,3*n)-0.5)*2*10^10; %Initial Positions for each body [m]
v0_all = (rand(1,3*n)-0.5)*2*10^4; %Initial Velocities for each body [m/s]

x0 = [r0_all, v0_all]'; %Initial States

xdot = dynamics_nbody(1, x0, const);

time = linspace(0, 3600*24*30, 2000); %Propagate over 1 month

[T, X] = ode45(@dynamics_nbody, time, x0, odeset('AbsTol',1e-11,'RelTol',1e-9), const);
