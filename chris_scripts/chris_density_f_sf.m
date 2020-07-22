function [eta, del] = chris_density_f_sf(rho)

rho0 = 1.784;
rho1 = 0.1786;

C = (rho-rho0)/(rho1-rho0);

[i5,~,i95] = find_5_50_95(C);

del = abs(i95-i5);

[nx, ny] = size(rho);

rm = imresize(mean(rho')',[nx ny]);

rhop = rho-rm;

[drpdx,~] = chris_gradient(rhop,1,1);

F = drpdx./rho.^2;

G = zeros(ny,nx);

for i = 1:ny
   G(i,:) = trapz(F(:,1:i)')'; 
end

eta = trapz(trapz(G)');

end

