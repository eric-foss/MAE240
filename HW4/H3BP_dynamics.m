function xdot = H3BP_dynamics(~, X)
%H3BP_DYNAMICS 

x = X(1);
y = X(2);
z = X(3);
xp = X(4);
yp = X(5);

r = sqrt(x^2 + y^2 + z^2);

xdot = zeros(6, 1);

xdot(1:3) = X(4:6);

xdot(4) = 2*yp - x/(r^3) + 3*x;
xdot(5) = -2*xp - y/(r^3);
xdot(6) = -z/(r^3) - z;

end

