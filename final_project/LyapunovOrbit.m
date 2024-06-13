
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%                              Main Function                              %
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X_LO, Phi_LO] = LyapunovOrbit(L, const)
%-------------------------------------------------------------------------%
%                          Function Description                           %
%-------------------------------------------------------------------------%

% The following user-defined function analyzes stable manifold trajectories
% about Lyapunov orbits within the PCR3BP. To do this, the function takes
% in a set of parameters defined in the structure array 'const' as well as
% the preferred colinear Lagrange point 'L1', 'L2', or 'L3' for which to do
% the analysis. Methods employed include eigenvector analysis, differential
% correction, and numerical integration techniques.



%-------------------------------------------------------------------------%
%                              Inputs Guide                               %
%-------------------------------------------------------------------------%

% L    :   a string containing the colinear Lagrange point to be used in
%          Lyapunov orbit determination; has the form 'L1', 'L2', or 'L3'.

% const:   a structure array containing the following values
%
% --> mu: the mass ratio of the system, this will be automatically
% calculated using the supplied primary masses.
%
% --> Rprimary: radius of the primary in kg
%
% --> Rsecondary: radius of the secondary in kg
%
% --> R: the distance between the primaries in km
%
% --> eps: the desired precision for the differential correction algorithm
%
% --> disp: the nondimension displacement from the desired Lagrange point
% to be used as the initial condition
%
% --> LOperiod: nondimensional guess for the period of the Lyaponov orbit
% (usually a Lyapunov orbit has a period half that of the primaries i.e.
% T_LO = pi)
%
% --> nm: number of stable manifold trajectories desired
%
% --> manprop: number of periods (of the primaries) to propagate the
% manifold trajectories backward
%
% --> I: max number of iterations for differential correction method



%-------------------------------------------------------------------------%
%                             Outputs Guide                               %
%-------------------------------------------------------------------------%

% Figure 1: a plot of the propigated initial guess without correction
% applied.

% Figure 2: a plot of the corrected periodic orbit with marked Lagrange
% point

% Figure 3: a plot of the normalized (J/J0) Jacobi constant vs time
% (non-dimensional)

% Figure 4: a plot of the stable and unstable manifold trajectories with the Lyapunov
% orbit and primaries




%-------------------------------------------------------------------------%
%         Symbolic Determination of Important Algorithm Equations         %
%-------------------------------------------------------------------------%

% Define symbolic variables
syms x y xdot ydot mu

% Define the effective potential for the PCR3BP
V = 1/2*(x^2 + y^2) + (1 - mu)/sqrt((x + mu)^2 + y^2) +...
    mu/sqrt((x - 1 + mu)^2 + y^2);

% Define the gradient
Vx = diff(V, x);
Vy = diff(V, y);

% Determine equations of motion
xddot = 2*ydot + Vx;
yddot = -2*xdot + Vy;

% Define state space form of the dynamics
F = [xdot; ydot; xddot; yddot];

% Define the state vector
X = [x; y; xdot; ydot];

% Determine the A matrix
A = simplify(jacobian(F, X));

% Convert helpful equations for use in functions
eqns.A = matlabFunction(A);
eqns.xddot = matlabFunction(simplify(xddot));
eqns.yddot = matlabFunction(simplify(yddot));



%-------------------------------------------------------------------------%
%                        Helpful Numerical Values                         %
%-------------------------------------------------------------------------%

% Options for use in ODE45 and fzero functions
const.optionsODE = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
const.optionsZERO = optimset('DiffMaxChange', const.eps, 'TolCon', ...
    const.eps, 'TolFun', const.eps);



%-------------------------------------------------------------------------%
%                  Colinear Lagrange Point Specification                  %
%-------------------------------------------------------------------------%

% Determine the polynomial associated with the lagrange point location
% given user specified input
switch L
    case 'L1'
        
        % Define the polynomial for L1
        p = [1, (3 - const.mu), (3 - 2*const.mu), const.mu, 2*const.mu, const.mu];

    case 'L2'
        
        % Define the polynomial for L2
        p = [1, (3 - const.mu), (3 - 2*const.mu), -const.mu, -2*const.mu, -const.mu];

    case 'L3'
        
        % Define the polynomial for L3
        p = [1, (3 - const.mu), (3 - 2*const.mu), (2 - const.mu), 2*const.mu, const.mu];

    otherwise

        % Return error if defined point does not fit the cases
        error("Error, invalid point specified. Use 'L1', 'L2', or 'L3'.")      
end

% Determine the real root of the polynomial
rts = roots(p);
for i = 1:5
    if isreal(rts(i))

        % Save the distance from the primary
        gamma = rts(i);

    end
end

% Detemine the x location of the defined lagrange point
xL = (1 - const.mu) + gamma;
XL = [xL; 0; 0; 0];



%-------------------------------------------------------------------------%
%                       Initial Condition Generation                      %
%-------------------------------------------------------------------------%

% Determine the eigenvalues and vectors associated with the linearized
% dynamics
[e_vec, e_val] = eigs(eqns.A(const.mu, xL, 0));

% Sort through eigenvalues to determine the imaginary one
LO_dir = zeros(4, 1);
flag = true;
for i = 1:length(e_val)

    % Determine if egienvector is imaginary
    if ~isreal(e_val(i, i)) && flag

        LO_dir = e_vec(:, i);
        
        % Sort through imaginary eigenvector
        for j = 1:length(e_val)
            
            % Pull out the real components
            LO_dir(j) = abs(real(LO_dir(j)));

        end
        
        % Ensure positive x displacement and negative y velocity
        LO_dir = LO_dir.*[1/LO_dir(1); 0; 0; -1/LO_dir(1)];

        % Set flag to false now that imaginary eigenvector has been found
        flag = false;
    end

end

% Determine the initial condition for the Lyapunov orbit
X0 = XL + const.disp*LO_dir;



%-------------------------------------------------------------------------%
%                   Pre-correction Orbit Determination                    %
%-------------------------------------------------------------------------%

% Integrate the equations of motion based on initial guess
tspan_guess = 0:0.01:const.LOperiod;
[~, X_guess] = ode45(@(t, X) PCR3BP_dyn(t, X, const, eqns), tspan_guess, X0, const.optionsODE);

% Plot the Pre-correction orbit
figure(1)
plot(X_guess(:, 1), X_guess(:, 2), 'k', X_guess(1, 1), X_guess(1, 2),'rv', xL, 0, 'b*',... 
    'MarkerSize', 10, 'LineWidth', 2), grid, xlabel('X/R'), ylabel('Y/R'), ...
    title('Pre-correction Trajectory')



%-------------------------------------------------------------------------%
%                   Post-correction Orbit Determination                   %
%-------------------------------------------------------------------------%

% Define the correction
[T, X_corrected] = LyapunovCorrect(X0, const, eqns);
tspan_corrected = 0:0.001:T;
if length(tspan_corrected) < 2
        error(['Error in Lyapunov orbit correction, no viable periodic orbit found. ' ...
            'Consider choosing a displacement closer to the desired Lagrange point or' ...
            ' increasing the the guess for the orbit''s period.'])
end

% Integrate new initial condition
[t_LO, X_LO] = ode45(@(t, X) PCR3BP_dyn(t, X, const, eqns), tspan_corrected, X_corrected, const.optionsODE);

% Define components for easier use
x_LO = X_LO(:, 1);
y_LO = X_LO(:, 2);
v_LO = X_LO(:, 3:4);

% Check to ensure an orbit has been found
if ~(xL < max(x_LO) && xL > min(x_LO))
        error(['Error in Lyapunov orbit correction, no viable periodic orbit found. ' ...
            'Consider choosing a displacement closer to the desired Lagrange point or' ...
            ' increasing the the guess for the orbit''s period.'])
end

% Determine the Jacobi Integral
J_LO = 1/2.*dot(v_LO, v_LO, 2) - 1/2.*(x_LO.^2 + y_LO.^2) - ...
    (1 - const.mu)./sqrt((x_LO + const.mu).^2 + y_LO.^2) - ...
    const.mu./sqrt((x_LO - 1 + const.mu).^2 + y_LO.^2);

% Plot the Lyapunov Orbit
figure(2)
plot(x_LO, y_LO, 'k', x_LO(1), y_LO(1),'rv', xL, 0, 'b*',... 
    'MarkerSize', 10, 'LineWidth', 2), grid, xlabel('X/R'), ylabel('Y/R'), ...
    legend({['T =', num2str(T), '; J =', num2str(J_LO(1))], '', L}), ...
    title(['Post-correction Trajectory (Lyapunov Orbit about ', L, ')'])

% Plot the Jacobi constant normalized about initial value
figure(3)
plot(t_LO, J_LO./J_LO(1), 'LineWidth', 2), grid, xlabel('tN'), ...
    ylabel('J/J_0'), title('Non-dimensional Jacobi Constant (Normalized)')



%-------------------------------------------------------------------------%
%                Stable Manifold Trajectory Determination                 %
%-------------------------------------------------------------------------%

% Determine the correct increment size for manifolds
inc = floor(size(X_LO, 1)/const.nm);


% Define the Monodromy matrix
Phi_T = reshape(X_LO(end, 5:20), 4, 4);

% Find the stable and unstable eigenvalues associated with Monodromy matrix
[eig_vec, eig_val] = eigs(Phi_T);
stable_ind = 0;
unstable_ind = 0;
stable_val = 1e15;
unstable_val = -1e15;
for j = 1:length(eig_val)
    if isreal(eig_val(j, j)) && eig_val(j, j) < stable_val
        stable_val = eig_val(j, j);
        stable_ind = j;
    end
    if isreal(eig_val(j, j)) && eig_val(j, j) > unstable_val
        unstable_val = eig_val(j, j);
        unstable_ind = j;
    end
end

% Ensure stable and unstable eigenvector exist
if stable_ind == 0 || unstable_ind == 0
    error(['Error in manifold trajectory generation, no stable eigenvectors exist.' ...
        ' Consider choosing a displacement closer to the desired Lagrange point.'])
end
% Define the stable and unstable eigenvectors
stable_vec = eig_vec(:, stable_ind)/norm(eig_vec(:, stable_ind));
unstable_vec = eig_vec(:, unstable_ind)/norm(eig_vec(:, unstable_ind));

% Determine manifold trajectory for each iteration
check = 0;
for s = -1:2:1
    for i = 1:const.nm
        
        % Determine currend index desired
        k = i*inc;
        
        % Apply a small purtubation along the stable eigenvector direction
        X0_stable = X_LO(k, 1:4)' + s*1e-5*stable_vec;
    
        % Apply a small purtubation along the stable eigenvector direction
        X0_unstable = X_LO(k, 1:4)' + s*1e-5*unstable_vec;
        
        % Determine backwards time span for stable integration
        tspan_stable = 0:-0.01:-const.manprop*T;
    
        % Determine forward time span for unstable integration
        tspan_unstable = 0:0.01:const.manprop*T;
    
        % Integrate to determine manifold trajectories
        [~, X_stable] = ode45(@(t, X) PCR3BP_dyn(t, X, const, eqns), tspan_stable, X0_stable, const.optionsODE);
        [~, X_unstable] = ode45(@(t, X) PCR3BP_dyn(t, X, const, eqns), tspan_unstable, X0_unstable, const.optionsODE);
        
        % Plot the stable and unstable manifold trajectories
        figure(4); hold on
        plot(X_stable(:, 1), X_stable(:, 2), 'b', X_unstable(:, 1), X_unstable(:, 2), 'r')

        if s == -1
            for j = 1:400
                
                dist = sqrt((X_unstable(j, 1) - (1 - const.mu))^2 + (X_unstable(j, 2))^2);
                
                if dist <= (const.Rsecondary+600)/const.R
                    figure(6); hold on
                    plot(X_unstable(1:400, 1), X_unstable(1:400, 2), 'r');
                    check = 1;
                    break;
                end

            end

            figure(5); hold on
            plot(X_unstable(1:400, 1), X_unstable(1:400, 2), 'r');

        end

    end
end


% Figure 4
figure(4); hold on

% Plot Saturn and Titan
viscircles([-const.mu 0], const.Rprimary/const.R, 'Color', '[0.8500 0.3250 0.0980]');
viscircles([(1 - const.mu) 0], const.Rsecondary/const.R, 'Color', '[0.6350 0.0780 0.1840]');

% Replot the Lyapunov orbit
plot(X_LO(:, 1), X_LO(:, 2), 'k', 'LineWidth', 2), grid on, xlabel('X/R'), ylabel('Y/R'), ...
    title('Stable and Unstable Manifold Trajectories')

legend({'Stable', 'Unstable'})
hold off
axis equal


% Figure 5
figure(5); hold on

% Plot Saturn and Titan
viscircles([-const.mu 0], const.Rprimary/const.R, 'Color', '[0.8500 0.3250 0.0980]');
viscircles([(1 - const.mu) 0], const.Rsecondary/const.R, 'Color', '[0.6350 0.0780 0.1840]');

% Replot the Lyapunov orbit
plot(X_LO(:, 1), X_LO(:, 2), 'k', 'LineWidth', 2), grid on, xlabel('X/R'), ylabel('Y/R'), ...
    title('Quick Descent Unstable Manifolds')

hold off
xlim([0.98 1.06]);
ylim([-0.03 0.03]);
daspect([1 1 1]);

% Figure 6
figure(6); hold on

% Plot Saturn and Titan
viscircles([-const.mu 0], const.Rprimary/const.R, 'Color', '[0.8500 0.3250 0.0980]');
viscircles([(1 - const.mu) 0], const.Rsecondary/const.R, 'Color', '[0.6350 0.0780 0.1840]');

% Replot the Lyapunov orbit
plot(X_LO(:, 1), X_LO(:, 2), 'k', 'LineWidth', 2), grid on, xlabel('X/R'), ylabel('Y/R'), ...
    title('Quick Descent Unstable Manifolds')

hold off
xlim([0.98 1.06]);
ylim([-0.03 0.03]);
daspect([1 1 1]);

if check == 0
    disp('No Quick Descent Unstable Manifolds Exist For this Orbit. Consider a farther displacement');
end


% Reshape state vector for output
Phi_LO = zeros(4, 4, length(t_LO));
for i = 1:length(t_LO)
    Phi_LO(:, :, i) = reshape(X_LO(i, 5:end), 4, 4);
end
X_LO = X_LO(:, 1:4);

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%                            Auxilary Functions                           %
%-------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------------------------------------------------------------------%
%                            Dynamics Function                            %
%-------------------------------------------------------------------------%

function Xdot = PCR3BP_dyn(~, X, const, eqns)
% The following user-defined function outlines the dynamics of the Planar
% Circular Restricted Three Body Problem for use in numerical integration.
% The equations of motion are in the non-dimensional form.

% Define state terms
x = X(1);
y = X(2);
xdot = X(3);
ydot = X(4);

% Determine accelerations from given eqns of motion
xddot = eqns.xddot(const.mu, x, y, ydot);
yddot = eqns.yddot(const.mu, x, xdot, y);

% Determine desired output based on user input
% If state greater than 4, the STM is assumed to be appended to the state
% vector
if length(X) > 4

    % Determine current STM
    Phi0 = reshape(X(5:end), 4, 4);

    % Find current A matrix
    A = eqns.A(const.mu, x, y);

    % Combine propgated state and STM
    Xdot = [xdot; ydot; xddot; yddot; reshape(A*Phi0, 16, 1)];

else

    % Otherwise, return only the progogated state
    Xdot = [xdot; ydot; xddot; yddot];
  
end
end



%-------------------------------------------------------------------------%
%                   Differential Correction Function                      %
%-------------------------------------------------------------------------%

function [T, X0_c] = LyapunovCorrect(X0, const, eqns)
% the following user-defined function uses differential correction in an
% iterative procedure to refine Lyupanov orbit intitial conditions. Such
% orbits are expected to have intitial conditions that lay within the XY
% plane of the synodic CR3BP frame and have trajectories perpendicular to
% the X-axis.

% Define loop variables
eps = const.eps; 
I = const.I; 
xdot_cross = const.eps + 1;
i = 0;

% Define the initial condition with STM
X0 = [X0; reshape(eye(4), 16, 1)];

% Define convergence loop to run until x velocity is below convergence
% criteria
while abs(xdot_cross) > const.eps && i < I

    % Update counter
    i = i + 1;

    % Determine where and when initial condition trajectory will recross
    % the x-axis
    [t_cross, X_cross] = xCrossing(X0, const, eqns);

    % Define the state variables at crossing
    x_cross = X_cross(1);
    y_cross = X_cross(2);
    xdot_cross = X_cross(3);
    ydot_cross = X_cross(4);
    
    % Define the STM at the time of crossing
    phi_cross = reshape(X_cross(5:end), 4, 4);

    % Determine the x acceleration at crossing for use in later calc
    xddot_cross = eqns.xddot(const.mu, x_cross, y_cross, ydot_cross);

    % Determine the deviation in the initial y velocity for this iteration
    % by holding the initial x position constant
    dydot0 = -xdot_cross/(phi_cross(3,4) - xddot_cross/ydot_cross*phi_cross(2,4));

    % Apply the correction to the initial state
	X0(4) = X0(4) + dydot0;

end

% Determine period of the orbit
T = 2*t_cross;

% Define new initial state for periodic orbit
X0_c = X0;

end



%-------------------------------------------------------------------------%
%                        X-axis Crossing Fucntion                         %
%-------------------------------------------------------------------------%

function [t_crossing, X_crossing] = xCrossing(X0, const, eqns)
% The following user-defined function determines where and when a trajectory
% will recross the axis of primaries within the CR3BP.

% Define global variable for passing between xCrossing and zeroFunc
% This is to avoid unnecessary reintigration of the equations of motion
% X_zero is the state at the crossing of the x axis
global X_zero

% Define initial guess to be slightly less than 1 half period of the
% standard Lyapunov orbit (roughly when the recrossing would occur)
% NOTE: 2*pi --> period of primaries orbit
% NOTE: pi   --> period of standard Lyapanov
t_before = const.LOperiod/2;

% Integrate the initial condition to determine the state just before
% crossing
[~, Xx] = ode45(@(t, X) PCR3BP_dyn(t, X, const, eqns), [0, t_before], X0, const.optionsODE);
X_before = Xx(end, :);

% Determine when and where crossing occurs
t_crossing = fzero(@(t) zeroFunc(t, t_before, X_before, const, eqns), t_before, const.optionsZERO);
X_crossing = X_zero;
clear global X_zero

end



%-------------------------------------------------------------------------%
%                            Find Zero Function                           %
%-------------------------------------------------------------------------%

function y_crossing = zeroFunc(t, t_before, X_before, const, eqns)
% The following user-defined function returns the y position of a crossing
% so that it may be used to converge onto the exact point at which the
% crossing occurs.
global X_zero

% Determine the state associated with x-axis crossing
if t == t_before

    X_zero = X_before;

else
   
    % Integrate the equations of motion to determine y pos at end state
    [~, Xs_zero] = ode45(@(t, X) PCR3BP_dyn(t, X, const, eqns), [t_before, t], X_before, const.optionsODE);
    X_zero = Xs_zero(end, :);

end

% Return y position for use in zero location determination
y_crossing = X_zero(2);

end