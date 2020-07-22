function [Xi,dcdx,dcdy] = chris_scalar_dissipation(C,dpxdx)

    D = 1.5*10^(-5);
    [dcdx,dcdy] = chris_gradient(C,1,1);
    Xi = D*(dpxdx.^2)*dcdx.^2+dcdy.^2; 

end

