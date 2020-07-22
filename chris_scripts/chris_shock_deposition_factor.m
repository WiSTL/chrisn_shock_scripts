function [eta] = chris_shock_deposition_factor(rho)

[~,drdy] = chris_gradient(rho,1,1);

eta = nansum(nansum(abs(medfilt2(drdy))./rho.^2));

end

