function [P,rho,T,P0,T0,M2] = moving_shock(M,g)

P = 1 + (2*g/(g+1))*(M.^2 - 1);

T = P.*(((g+1)/(g-1))+P)./(1+((g+1)/(g-1))*P);

rho = P./T;

M2 = (2/(g+1))*(M+1./M).*sqrt(1./T);

P0 = P.*(1+((g-1)/2)*M2.^2).^(g/(g-1));

T0 = T.*(1+((g-1)/2)*M2.^2);

end

