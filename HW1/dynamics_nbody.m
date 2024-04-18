function xdot = dynamics_nbody(~, x, const)
%DYNAMICS_NBODY

n = length(x)/6 - 1; %number of bodies

r = zeros(3, n); %
a = zeros(3, n);

for i = 1:n
    r(:, i) = x(3*i-2:3*i);
end

R = x(length(x)/2 - 2:length(x)/2);

for i = 1:n
    for j = 1:n
        if i == j
            continue;
        end
        rji = r(:, i) - r(:, j);
        a(:, i) = a(:, i) - const.G*const.m_all(j)*(rji)/norm(rji)^3;
    end
end

A = zeros(3, 1);

for i = 1:n
    rjN1 = R - r(:, i);
    A = A - const.G*const.m_all(i)*rjN1/norm(rjN1)^3;
end

acc = [reshape(a, 1, 3*n)'; A];

xdot = [x(length(x)/2 + 1:end); acc];