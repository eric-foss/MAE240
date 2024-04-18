function xdot = dynamics_nbody_bodycentric(~, x, const)
%DYNAMICS_NBODY_BODYCENTRIC

n = height(x)/6;

a = zeros(3, n-1);



r = reshape(x(1:3*(n-1)), 3, n-1);
R = x(3*n - 2: 3*n);

for i = 1:(n-1)   
    a(:, i) = -const.G*(const.m_all(n)+const.m_all(i))*r(:, i)/norm(r(:, i))^3;
    for j = 1:(n-1)
        if j == i
            continue;
        end
        rji = r(:, i) - r(:, j);
        a(:, i) = a(:, i) - const.G*const.m_all(j)*(rji/(norm(rji)^3) + r(:, j)/(norm(r(:, j))^3));
    end
end


A = -const.G*const.m_all(n)*R/norm(R)^3;
for j = 1:(n-1)
    rjN1 = R - r(:, j);
    A = A - const.G*const.m_all(j)*(rjN1/norm(rjN1)^3 + r(:, j)/norm(r(:, j))^3);
end

acc = [reshape(a, 3*(n-1), 1); A];

xdot = [x(3*n+1: end); acc];

end
