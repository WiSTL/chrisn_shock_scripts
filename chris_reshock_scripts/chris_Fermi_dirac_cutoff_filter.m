function [um,vm,rhom,up,vp,rhop] = chris_Fermi_dirac_cutoff_filter(u,v,rho,C)

[l_ux,l_uy,l_vx,l_vy,l_cx,l_cy] = chris_taylor_microscale(u,v,C);

l_u = sqrt(l_ux^2+l_uy^2);
l_v = sqrt(l_vx^2+l_vy^2);
l_c = sqrt(l_cx^2+l_cy^2);

lambda_t = mean([l_u,l_v,l_c])

[uf,kx,kz] = fft2d(u);
[vf] = fft2d(v);
[rhof] = fft2d(rho);

lambda_l = 0.75*lambda_t;
kc = 2*pi/lambda_l;
[Kx,Kz] = meshgrid(kx,kz);
Kxz = sqrt(Kx.^2 + Kz.^2);

F = (1+exp((Kxz-kc)/(0.1*kc))).^(-1);

umf = uf.*F;
vmf = vf.*F;
rhomf = rhof.*F;

um = ifft2d(umf);
vm = ifft2d(vmf);
rhom = ifft2d(rhomf);

up = u-um;
vp = v-vm;
rhop = rho-rhom;
end

