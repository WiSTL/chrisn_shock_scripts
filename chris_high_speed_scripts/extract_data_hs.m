function [out] = extract_data_hs(C,params)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% run 11 params 
% dt = 5E-5
% hd0 = 28.494
% h0 = 0.00512
% ts = 0.00025
% r1 = 5.096
% r2 = 7.957
% nu = 7.011914E-06
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

out = {};

% C11 = C_corrected(500:1080,160:502,5:95)/0.7;
C(C<0)=0;
C(C>1)=1;

out.C = C;

[d,ds] = chris_exp_length_scales_hs(C);
%[xl,A] = chris_exp_pdfs_hs(C11(:,:,7:65));

[nz,nx,nt] = size(C);

dt = params(1);

D = 0.00059; %m2/s wasik and McCulloh 1968

dpxdx = 4081.63;

hd0 = params(2);
h0 = params(3);
out.h0 = h0;

t0 = h0/hd0;
ts = params(4);

out.h0 = h0;
out.ts = ts;
out.t0 = t0;

tcl = t0*(100/5)*(1.6*10^5)^(-0.5);
tcu = t0*(100/sqrt(15))*(1.6*10^5)^(-0.5);

r1 = params(5);
r2 = params(6);

out.rho_1 = r1;
out.rho_2 = r2;

nu = params(7);

out.nu = nu;

out.rho = r1 + (r1-r2)*C;

out.MM = out.rho.*C.*(1-C);

out.rho_m = mean(out.rho,2);
out.C_m = mean(C,2);

MM11_norm = squeeze(out.rho_m.*out.C_m.*(1-out.C_m));

psi = sum(squeeze(sum(out.MM,1)),1)./(sum(MM11_norm,1)*nx);


for n = 1:nt
    Cm(:,n) = mean(squeeze(C(:,:,n)),2);
    
    [~,out.eta(:,n),x(:,n),out.delta(n),out.y0(n)] = chris_C(squeeze(C(:,:,n)));
    
    for i = 1:nx
    out.fluc.Cp(:,i,n) = C(:,i,n) - Cm(:,n);
    end
    
    out.fluc.cp2(:,n) = mean(squeeze(out.fluc.Cp(:,:,n).^2),2);
    out.fluc.cp3(:,n) = mean(squeeze(out.fluc.Cp(:,:,n).^3),2);
    out.fluc.cp4(:,n) = mean(squeeze(out.fluc.Cp(:,:,n).^4),2);
end

x = [1:nx]/dpxdx;
out.t = [0:nt-1]*dt;
out.tau = (out.t-out.ts)/out.t0;
% e = eta(:,1);
% [E,T] = meshgrid(e,t);

[E,Z,X,T] = meshgrid_hs(out.eta,x,out.t,out.delta);

out.grid.E = E;
out.grid.Z = Z;
out.grid.X = X;
out.grid.T = T;

h = out.delta/dpxdx;
out.hd = chris_derivative(smooth(h),mean(diff(out.t)));
out.hdd = chris_derivative(out.hd,mean(diff(out.t)));

out.V0 = chris_derivative(out.y0/dpxdx,mean(diff(out.t)));
out.a0 = chris_derivative(out.V0,mean(diff(out.t)));

out.Re_h = (h).*out.hd'/nu;
out.C_h = (out.hd./(h)')*t0;

[out.dcz,out.dcx,out.dct] = chris_gradient_hs(C,X,Z/dpxdx,T);
[out.dcpz,out.dcpx,out.dcpt] = chris_gradient_hs(out.fluc.Cp,X,Z/dpxdx,T);

out.fluc.xi = 2*D*squeeze(mean(out.dcpz.^2+out.dcpx.^2,2));
out.fluc.xi_x = squeeze(mean(out.dcpx.^2,2));
out.fluc.xi_z = squeeze(mean(out.dcpz.^2,2));

lambda_x = sqrt(out.fluc.cp2./out.fluc.xi_x);
lambda_z = sqrt(out.fluc.cp2./out.fluc.xi_z);

lambda_z(isnan(lambda_z))=0;
lambda_z(isinf(lambda_z))=0;
lambda_x(isnan(lambda_x))=0;
lambda_x(isinf(lambda_x))=0;

lambda = sqrt(lambda_x.^2 + lambda_z.^2);

lambda_x_m = trapz([1:nz]/dpxdx,lambda_x)./h;
lambda_z_m = trapz([1:nz]/dpxdx,lambda_z)./h;

out.fluc.xi_m = trapz([1:nz]/dpxdx,out.fluc.xi)./h;
out.fluc.cp2_m = trapz([1:nz]/dpxdx,out.fluc.cp2)./h;

out.fluc.C_xi = (h./out.hd').*out.fluc.xi_m./out.fluc.cp2_m;

for i = 1:nt
    out.fluc.dwc(:,i) = -1*squeeze(mean(out.dct(:,:,i),2))+out.V0(i)*squeeze(mean(out.dcz(:,:,i),2));
    out.fluc.wc(:,i) = cumtrapz(out.eta(:,i),out.fluc.dwc(:,i));
end

for i = 1:nx
    for j = 1:nt
        out.length_scales.h_spatial(i,j) = trapz(squeeze(C(:,i,j)).*(1-squeeze(C(:,i,j))))/dpxdx;
    end
end

out.length_scales.h = h;
out.length_scales.lambda_T_x = lambda_x;
out.length_scales.lambda_T_z = lambda_z;
out.length_scales.lambda_T_x_m = lambda_x_m;
out.length_scales.lambda_T_z_m = lambda_z_m;



end
