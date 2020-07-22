function [dudx] = chris_compact_derivative_8(u,dx)

[dudx] = chris_derivative(u,dx)';

nu = length(u)
f = [];
g = [];

%tridiag coefficients
a = ones(nu-5);%ones(n-4)*(-3/20);
b = ones(nu-5)*(-3/20);
c = b;

% interpolation
for i=3:nu-3
    f(i-2) = (-1/25)*(u(i+3)+u(i-2))-(61/100)*(u(i+2)+u(i-1)) + (u(i+1)+u(i));
end

f_h = tridiag( a, b, c, f );

n = length(f_h)

for i = 3:n-1
    g(i-2) = (2/dx)*(f_h(i) - f_h(i-1)) - (61/100)*(1/dx)*(u(i+1) - u(i-1)) - (2/75)*(1/dx)*(f_h(i+1) - f_h(i-2));
end

n = length(g);
a = ones(n-5);%ones(n-4)*(-3/20);
b = ones(n-5)*(-3/20);
c = b;

dudx1 = tridiag( a, b, c, g );

dudx(9:nu-4) = dudx1(3:n-2);

end