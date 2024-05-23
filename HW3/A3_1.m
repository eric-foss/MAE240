%% MAE 240 Assignment 3 Problem 1
% Eric Foss
% A17068006
clear; close all; clc;

%% SPICE SETUP

%  Add specified folder to the top of the search path
addpath c:\Users\Eric\MAE240\HW3\mice\mice\lib
addpath c:\Users\Eric\MAE240\HW3\mice\mice\src\mice

%  Load SPICE kernel files into MATLAB
cspice_furnsh('de430.bsp');                       % Planetary ephemerides
cspice_furnsh('pck00010.tpc');                    % Planet orientation and radii
cspice_furnsh('gm_de431.tpc');                    % GM/mass 
cspice_furnsh('naif0011.tls'); 

%% CONSTANTS/PARAMETERS

% Numerical Integration Paramters 
utc0 = '2020/11/04 00:00:00 UTC';  
const.et0 = cspice_str2et( utc0 );
jyear = cspice_jyear;
Nt = 30000;
time = linspace(const.et0, const.et0 + jyear, Nt);

% Gravitational parameters [km^3/s^2]
const.earth.mu = cspice_bodvrd( 'Earth' , 'GM' , 1 );  
const.moon.mu = cspice_bodvrd( 'Moon' , 'GM' , 1 );
const.nu = 0.0122;

% Dimensionless Scales
const.R = 384400; % Distance between earth and moon [km] (length scale)
const.N = sqrt((const.earth.mu + const.moon.mu)/const.R^3); % 1/Time scale [s]


%  Earth radii (mean equatorial, km)
R = cspice_bodvrd( 'Earth', 'RADII', 3);
const.earth.R = ( R(1) + R(2) ) / 2;


%% NUMERICAL INTEGRATION

r0 = [40000; 7000; 3000]/const.R;


v0 = [0; 3; 0]/(const.R*const.N);

x0 = [r0; v0];
tspan = time*const.N;

[t, X] = ode45(@CR3BP_dynamics, tspan , x0, odeset('AbsTol',1e-13,'RelTol',1e-13), const);

r = X(:, 1:3);
x = r(:, 1);
y = r(:, 2);
z = r(:, 3);
v = X(:, 4:6);


%% JACOBIAN

J = 0.5*dot(v, v) - 0.5*(x.^2 + y.^2) - (1 - const.nu)./sqrt((x + const.nu).^2 + y.^2 + z.^2) - const.nu./sqrt((x-1+const.nu).^2 + y.^2 + z.^2);

J_norm = zeros(length(J), 1);
J0 = norm(J(1, :));

for i = 1:length(J)
    J_norm(i) = norm(J(i, :));
end


%% 

[X_m, ~] = cspice_spkezr('Moon', time, 'ECLIPJ2000', 'none', 'Earth');
X_m = X_m';

X_inertial = zeros(size(X));

const.primary.mu = const.earth.mu;
const.secondary.mu = const.moon.mu;

for k = 1:length(t)

    moon = X_m(k, :)';
    sat = X(k, :)';

    X_inertial(k, :) = AGI_Qmethod_Inverse_RotateBarycenter_rv_to_J2000_rv(moon, sat, const)';

end

%% FIGURES

figure(1);
plot3(x, y, z);
xlabel('Dimensionless X'); ylabel('Dimensionless Y'); zlabel('Dimensionless Z');
title('Satellite Trajectory in Barycentric Synodic Frame');

figure(2);
plot3(X_inertial(:, 1), X_inertial(:, 2), X_inertial(:, 3));
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Satellite Trajectory in Ecliptic J2000 Inertial Frame')

figure(3);
plot(J_norm/J0);
xlabel('Dimensionless Time'); ylabel('Jacobian Magnitude');
title('Normalized Jacobian');