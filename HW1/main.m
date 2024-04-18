%% MAE 240 Assignment 1
% Eric Foss
% A17068006
clear; close all; clc;

%% Universal Constants
const.G = 6.6742e-11; %Gravitation Constant [N m^2/kg^2]

%% PROBLEM 1

n = input("How many bodies?");

const.m_all = rand(1, n+1) * 10^29; %random mass for each body [kg]
const.m_all(n+1) = const.m_all(n+1)/10^26; %(n+1)th mass negligible (but real)

r0_all = (rand(1,3*(n+1))-0.5)*2*10^10; %Initial Positions for each body [m]
v0_all = (rand(1,3*(n+1))-0.5)*2*10^4; %Initial Velocities for each body [m/s]

x0 = [r0_all, v0_all]'; %Initial States

time = linspace(0, 3600*24*30, 2000); %Propagated over 1 month

% NUMERICAL INTEGRATION (SEE dynamics_nbody.m)
[t, X] = ode45(@dynamics_nbody, time, x0, odeset('AbsTol',1e-11,'RelTol',1e-9), const);

% CONSERVATION OF ENERGY AND MOMENTUM FOR N BODIES

%Preallocating vectors
T = zeros(height(t), 1);
U = zeros(height(t), 1);
E = zeros(height(t), 1);
H = zeros(3, height(t));
h = zeros(height(t), 1);

%Looping through each stored state
for i = 1:height(X)
    
    %Position and Velocity matrix for current state
    for j = 1:n
        r(:, j) = X(i, 3*j-2:3*j);
        v(:, j) = X(i, 3*(n+1) + 3*j - 2:3*(n+1) + 3*j);
    end

    T(i) = 0.5*trace(v'*v.*const.m_all(1:n)); %Kinetic Energy for current state
    
    %Potential Energy/Angular Momentum Loop
    u = zeros(n, 1);
    for ii = 1:n
        for jj = 1:n
            if jj == ii
                continue;
            end
            
            u(ii) = u(ii) + const.G*const.m_all(ii)*const.m_all(jj)/norm(r(:, ii)-r(:, jj));
            
        end

        H(:, i) = H(:, i) + const.m_all(ii)*cross(r(:, ii), v(:, ii));

    end

    U(i) = 0.5*sum(u);  %Total Potential Energy

    E(i) = T(i) - U(i); %Total Energy

    h(i) = norm(H(:, i)); %Angular Momentum Magnitude

end

% FIGURES P1A

%Conservation of Energy
figure(1); hold on;
plot(t, E/E(1));
title('Conservation of Energy');
ylabel('Normalized Total Energy'); xlabel('Time (s)');

%Conservation of Momentum
figure(2); hold on;
plot(t, h/h(1));
title('Conservation of Momentum');
ylabel('Normalized Magnitude of Angular Momentum'); xlabel('Time (s)');


% CONSERVATION OF ENERGY AND MOMENTUM FOR (N+1)ST BODY

E_N1 = zeros(height(X), 1);
h_N1 = zeros(height(X), 1);

for i = 1:height(X)
    
    r = X(i, 3*n + 1: 3*(n+1));
    v = X(i, end-2:end);

    T = 0.5*const.m_all(end)*dot(v, v);

    U = 0;
    for j = 1:n
        
        rj = X(i, 3*j-2:3*j);
        rjN1 = r - rj;
        U = U + const.G*const.m_all(j)*const.m_all(end)/norm(rjN1);

    end

    E_N1(i) = T - U;

    H = const.m_all(end)*cross(r, v);
    h_N1(i) = norm(H);

end

% FIGURES P1B
%Conservation of Energy
figure(3); hold on;
plot(t, E_N1/E_N1(1));
title('Nonconservation of Energy of N+1st Body')
ylabel('Normalized Total Energy'); xlabel('Time (s)');

figure(4);
plot(t, h_N1/h_N1(1));
title('Nonconservation of Momentum of N+1st Body');
ylabel('Normalized Magnitude of Angular Momentum'); xlabel('Time (s)');


%% PROBLEM 2

rs0 = reshape(r0_all, 3, n+1);
rs0 = rs0 - rs0(:, n);
rs0(:, n) = [];

vs0 = reshape(v0_all, 3, n+1);
vs0 = vs0 - vs0(:, n);
vs0(:, n) = [];

xs0 = [reshape(rs0, n*3, 1); reshape(vs0, 3*n, 1)];

[ts, Xs] = ode45(@dynamics_nbody_bodycentric, time, xs0, odeset('AbsTol',1e-11,'RelTol',1e-9), const);

% COMPARISON TO P1
rs_f = reshape(Xs(end, 1:3*(n-1)), 3, n-1);
r_f = reshape(X(end, 1:3*n), 3, n);

rs_f_converted = rs_f + r_f(:, n);

rsN1_f = Xs(end, 3*n-2:3*n);
rN1_f = X(end, 3*(n+1)-2:3*(n+1));
rsN1_f_converted = rsN1_f + r_f(:, n)';