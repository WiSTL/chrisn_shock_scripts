function [a,rho,T,Mt,Mr] = chris_he_ar_interface_M16(C)

a = 436.9*exp(0.3089*C)+3.418*exp(4.387*C);
rho = 9.896*exp(0.04791*C)-6.44*exp(0.3974*C);
T = 565.3*exp(-0.09564*C)-0.2014*exp(4.668*C);
Mt = 1.835*exp(-0.01543*C)-0.04212*exp(0.8699*C);
Mr = 0.9622*exp(0.05632*C)+0.03772*exp(1.387*C);

end

