function [d0_dot] = chris_initial_growth(Ms,e,d0,k,g,a0,r1,r2)

alpha = (2*a0/(g+1))*(Ms^2 - 1)/Ms;

eta = exp(k*d0)-(log(r1)/log(r2));

betap = alpha*e*log(r2);
beta = betap*eta;

a = r1/r2;

c = (r2-r1)/(exp(1+e)*log(a));

d0_dot = 2*betap*exp(1+e)*((exp(k*d0)*(a)^2 - 1)/(2*log(a)+k*d0) + (1-a^2 + 2*(a^2)*log(a))*((log(r1)/log(r2))-1)/((log(a))^2) - ((a)^2 - 1)/(2*log(a)+k*d0))...
    -betap*r2*((exp(k*d0)*(a) - 1)/(log(a)+k*d0) + (1-a + 2*(a)*log(a))*((log(r1)/log(r2))-1)/((log(a))^2) - ((a) - 1)/(log(a)+k*d0))...
    - 2*beta*exp(1+e)*((a^2)*(k*d0*coth(k*d0)-2*log(a)) - (k*d0)/sinh(k*d0))/((k*d0)^2 - 4*(log(a))^2)...
    + beta*r2*((a)*(k*d0*coth(k*d0)-2*log(a)) - (k*d0)/sinh(k*d0))/((k*d0)^2 - 4*(log(a))^2);

d0_dot = d0_dot/c;

end

