clearvars; close all; %clc

q = [ (1/3)^(1/3) 0 0 ]; %given L2 equilibrium point
Uqq = [ (3-(1/norm(q)^3)+3*q(1)/norm(q)^4), (3*q(1)*q(2)/norm(q)^3),(3*q(1)*q(3)/norm(q)^5); ...
        (3*q(1)*q(2)/norm(q)^3), (-1/norm(q)^3 + 3*q(2)^2/norm(q)^5), (3*q(2)*q(3)/norm(q)^5); ...
        (3*q(1)*q(3)/norm(q)^5), (3*q(2)*q(3)/norm(q)^5), (-1-1/norm(q)^3 + 3*q(3)^2/norm(q)^5)]; % Uqq* matrix derived in lecture
z_hat = [ 0 -1 0; 1 0 0; 0 0 0 ]; %z_hat_tilde
zz_hat= -[ 1 0 0; 0 1 0; 0 0 0 ]; %z_hat_tilde dotted with z_hat_tilde in recording
U_bar = eye(3); %identity matrix
J = [zeros(3,3) eye(3); -eye(3) zeros(3,3)]; %jacobian?
Hxx = [(-zz_hat-Uqq) z_hat; -z_hat U_bar];
JHxx = J*Hxx;
sigma= eig(JHxx) %eigenvalues
[V,D]=eig(JHxx) %diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors, so that A*V = V*D


% %{
q = [(1/3)^(1/3); 0; 0]; 
prop_stable = V(4:6,1)*0.1; 
prop_unstable = V(4:6,2)*0.1;

% Initial conditions and orbit period
x01 = [vertcat(q, prop_stable)]';       % Stable Manifold (positive)
x02 = [vertcat(q, -prop_stable)]';      % Stable Manifold (negative)
x03 = [vertcat(q, prop_unstable)]';     % Unstable Manifold (positive)
x04 = [vertcat(q, -prop_unstable)]';    % Unstable Manifold (negative)

T = 2*pi*10.0;

% Integration time span
tvec = linspace(0, T, 5000);

% Set integrator options
tol = 1e-12;
options= odeset( 'RelTol' , tol , 'AbsTol' , tol );

% Integrate the equations of motion using ode45
[t1, x1] = ode45(@dynamics_HR3BP,flip(tvec), x01, options);
[t2, x2] = ode45(@dynamics_HR3BP,flip(tvec), x02, options); 
[t3, x3] = ode45(@dynamics_HR3BP,tvec, x03, options);
[t4, x4] = ode45(@dynamics_HR3BP,tvec, x04, options);

% 2D Plot
figure(1); hold all
    plot( x1(:,1) , x1(:,2) )
    plot( x2(:,1) , x2(:,2) )
    plot( x3(:,1) , x3(:,2) )
    plot( x4(:,1) , x4(:,2) )
    
    plot( x1(1,1) , x1(1,2) , 'k.' , 'MarkerSize' , 12.0 ) % IC
    plot( (1/3)^(1/3) , 0 , 'r+' , 'MarkerSize' , 12.0 )
    
%     xlim( [-1.5 , 1.5] )
%}   
  
    
