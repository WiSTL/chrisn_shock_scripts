function [rho,mu] = chris_rho(C)

    r1 = 4.6081;
    r2 = 7.622;    
    
    mu1 = 3.5*10^(-5);
    mu2 = 4.5*10^(-5);
    
    rho = r2+(r1-r2)*C;
    mu = mu2+(mu1-mu2)*C;

end

