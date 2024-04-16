%% MAE 240 Assignment 1 Problem 1
% Eric Foss
% A17068006
clear; close all; clc;

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

[t, X] = ode45(@dynamics_nbody, time, x0, odeset('AbsTol',1e-11,'RelTol',1e-9), const);

% CONSERVATION OF ENERGY AND MOMENTUM

T = zeros(height(t), 1);
U = zeros(height(t), 1);
E = zeros(height(t), 1);
H = zeros(3, height(t));
h = zeros(height(t), 1);

for i = 1:height(X)
    
    for j = 1:n
        r(:, j) = X(i, 3*j-2:3*j);
        v(:, j) = X(i, 3*n + 3*j - 2:3*n + 3*j);
    end

    T(i) = 0.5*trace(v'*v.*const.m_all);
    
    u = zeros(5, 1);
    for ii = 1:n
        for jj = 1:n
            if jj == ii
                continue;
            end
            
            u(ii) = u(ii) + const.G*const.m_all(ii)*const.m_all(jj)/norm(r(:, ii)-r(:, jj));
            
        end

        H(:, i) = H(:, i) + const.m_all(ii)*cross(r(:, ii), v(:, ii));

    end

    U(i) = 0.5*sum(u);

    E(i) = T(i) - U(i);

    h(i) = norm(H(:, i));

end


figure(1);
plot(t, E);

figure(2);
plot(t, h);




