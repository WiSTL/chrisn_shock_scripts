function [v,rho,p,dpdx,T,x,delta] = chris_shock_structure(M1)

x = -20:0.001:20;

a = 0.6356 + 0.3632*M1.^(-2.084);
b = 0.3187 - 0.4198*M1.^(-3.453);
c = -1.742 + 1.719*M1.^0.6487;

v = a - b*erf(c*x);
p = (1./v).*(1+((M1^2)/3).*(1-v.^2));
dpdx = diff(p)./diff(x);
rho = (1./v);
T = 1+(M1.^2).*(1-v.^2);

x1 = find(dpdx>0.005*max(dpdx));

delta = x(x1(end))-x(x1(1));

end

