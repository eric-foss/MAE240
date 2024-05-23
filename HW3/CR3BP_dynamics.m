function xdot = CR3BP_dynamics(~, X, const)
%CR3BP_DYNAMICS 

nu = const.nu;
x = X(1);
y = X(2);
z = X(3);
xp = X(4);
yp = X(5);


xdot = zeros(6, 1);

xdot(1:3) = X(4:6);

xdot(4) = 2*yp + x - ((1-nu)*(x+nu)/((x+nu)^2 + y^2 + z^2)^(3/2)) - (nu*(x-1+nu)/((x-1+nu)^2 + y^2 + z^2)^(3/2));

xdot(5) = -2*xp + y*(1 - (1-nu)/((x+nu)^2 + y^2 + z^2)^(3/2) - nu/((x-1+nu)^2 + y^2 + z^2)^(3/2));

xdot(6) = z*((nu-1)/((x+nu)^2 + y^2 + z^2)^(3/2) - nu/((x-1+nu)^2 + y^2 + z^2)^(3/2));







end

