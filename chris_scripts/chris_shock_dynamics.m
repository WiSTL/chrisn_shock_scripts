function [M] = chris_shock_dynamics(C,MWa,MWb,ga,gb,M0)

% ga = 1.503;
% gb = 1.664;
% 
% MWa = 8.059;
% MWb = 39.95;

dMW = MWa - MWb;

xi = C;

g = ((((1-ga)/(1-gb))*ga-gb)*xi + gb)./((((1-ga)/(1-gb)) - 1)*xi + 1);

dgdxi = (((1-ga)/(1-gb))*(ga-gb))./((((1-ga)/(1-gb))-1)*xi+1).^2;

MW = dMW*xi + MWb;

n = length(xi);

M = zeros(1,n);

M(1) = M0;

for i = 1:n-1
    
    dxi = xi(i+1) - xi(i);
    
    mu(i) = sqrt(((g(i)-1)*(M(i)^2)+2)/(2*g(i)*(M(i)^2)-(g(i)-1)));
    
    gf(i) = 1 + (2*mu(i)*(M(i)^2 + 2))/((g(i)-1)*(M(i)^2)+2);
    
    ff(i) = (gf(i)/(g(i)*(g(i)+1)))*(mu(i)-g(i));
    
    lambda(i) = (1+(2/(g(i)+1))*((1-mu(i)^2)/mu(i)))*(1+2*mu(i)+1/M(i)^2);
    
    dM(i) = ((1-M(i)^2)/(lambda(i)*M(i)^2))*((ff(i) + 0.5*gf(i)/g(i))*dgdxi(i) - dMW*gf(i)/MW(i)^2)*dxi/2;
    
    M(i) = (M(i) + dM(i));
    
    mu(i) = sqrt(((g(i)-1)*(M(i)^2)+2)/(2*g(i)*(M(i)^2)-(g(i)-1)));
    
    gf(i) = 1 + (2*mu(i)*(M(i)^2 + 2))/((g(i)-1)*(M(i)^2)+2);
    
    ff(i) = (gf(i)/(g(i)*(g(i)+1)))*(mu(i)-g(i));
    
    lambda(i) = (1+(2/(g(i)+1))*((1-mu(i)^2)/mu(i)))*(1+2*mu(i)+1/M(i)^2);
    
    dM(i) = ((1-M(i)^2)/(lambda(i)*M(i)^2))*((ff(i) + 0.5*gf(i)/g(i))*dgdxi(i) - dMW*gf(i)/MW(i)^2)*dxi/2;
    
    M(i+1) = M(i) + dM(i);
    
end