function [b] = chris_density_self_correlation(rho)

    [ny,nx] = size(rho);
    rho_p=rho;
    rho_div_p=rho;
    
    rho_m = mean(rho');
    rho_div_m = mean((1./rho)');
    
    for j = 1:nx
        rho_p(:,j) = rho(:,j) - rho_m';
    end
    
    for j = 1:nx
        rho_div_p(:,j) = 1./rho(:,j) - rho_div_m';
    end
    
    b = -1*mean((rho_div_p.*rho_p)');


end

