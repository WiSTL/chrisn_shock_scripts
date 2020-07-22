
eta = -3:0.001:3;

n = length(eta);

R = 0.6;
a = sqrt(exp(1)/2);
re = 10;

rm = 1+(R-1)*(0.5*erfc(eta));

am = 1./rm;

dvp2 = 2*eta.*exp(-2*eta.^2).*((a^2)*(1-2*eta.^2) - (2*a*am/re).*(2*eta.^2 - 3));

de = abs(eta(2)-eta(1));

vp2 = zeros(1,n);

for i=2:n
    vp2(i) = vp2(i-1) + de*dvp2(i);
end
