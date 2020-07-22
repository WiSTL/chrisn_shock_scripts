function [P,rho,P0,T,M2] = normal_shock(M,g)

M2 = sqrt((1+((g-1)/2).*M.^2)./(g.*M.^2 - (g-1)/2));

rho = ((g+1).*M.^2)./(2+(g-1).*M.^2);

P = 1 + (2*g./(g+1)).*(M.^2 - 1);

T = P/rho;

%P0 = ((((g+1)*M.^2)./(2+(g+1)*M.^2)).^(g/(g-1))).*(((g+1)./(2*g*M.^2 - (g-1))).^(1/(g-1)));

P0 = ((1+g.*M.^2)./(1+g.*M2.^2)).*((1+((g-1)/2).*M2.^2)./(1+((g-1)/2).*M.^2)).^(g./(g-1));

end

