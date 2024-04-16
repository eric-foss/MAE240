function xdot = dynamics_nbody(~, x, const)
%DYNAMICS_NBODY Summary of this function goes here
%   Detailed explanation goes here

n = length(x)/6; %number of bodies

xdot = zeros(length(x), 1);

r = zeros(3, n); %
a = zeros(3, n);
v = zeros(3, n);

for i = 1:n
    r(:, i) = x(3*i-2:3*i);
end

for i = 1:n
    for j = 1:n
        if i == j
            continue;
        end
        rji = r(:, j) - r(:, i);
        a(:, i) = a(:, i) - const.G*const.m_all(j)*(rji)/norm(rji)^3;
    end
end



xdot = [x(length(x)/2 + 1:end); reshape(a, 1, 3*n)'];