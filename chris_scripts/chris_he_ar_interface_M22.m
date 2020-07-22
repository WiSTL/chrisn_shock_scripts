function [a,rho,T,Mt,Mr] = chris_he_ar_interface_M22(C)

a = 572.7*exp(0.2556*C)+4.543*exp(4.249*C);
rho = 6.743*exp(-0.05773*C)-2.105*exp(0.888*C);
T = 972.5*exp(-0.19*C)-2.041*exp(3.625*C);
Mt = 2.708*exp(-0.02399*C)-0.09*exp(0.8551*C);
Mr = 0.9765*exp(0.102*C)+0.02329*exp(2.022*C);

end

