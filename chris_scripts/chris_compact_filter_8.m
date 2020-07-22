function [f_h] = chris_compact_filter_8(f)

f_h = f;

al = (-3/20);

a(1) = (93+al*70)/128;
a(2) = (7+al*18)/16;
a(3) = (14*al - 7)/32;
a(4) = (1/16) - al/8;
a(5) = al/64 - 1/128;


n = length(f);

g = zeros(n-10);

for i = 5:n-5
    for j = 1:5
        g(i-4) = g(i-4) + a(j)*(f(i-(j-1))+f(i+(j-1)));
    end
end

ni = length(g);

a = ones(ni-5);%ones(n-4)*(-3/20);
b = ones(ni-5)*(al);
c = b;


ff = 0.5*tridiag( a, b, c, g );
f_h(5:n-6) = ff;

end