function [psi,ncount,err,a,b,gamma] = chris_velocity_hs_4(C,x,eta,tau,h,L)

eta_i = -2:0.05:2;

[C_i] = chris_interp_hs(C,eta,tau,eta_i,tau);

C_i(isnan(C_i))=0;

lnh=log(h/h(1));

[ne,nx,nt] = size(C_i);

[~,etam] = meshgrid(tau,eta_i);

[E,~,X,T] = meshgrid_hs(etam,x/L,lnh,h);

[dce,dcx,dct] = chris_gradient_hs(C_i,X,E,T);

dce = medfilt3(dce);
dcx = medfilt3(dcx);
dct = medfilt3(dct);

for k = 1:nt
    for j = 1:ne
        gamma(j,:,k) = (eta_i(j).*dce(j,:,k) - dct(j,:,k)).*L./h(k);
    end
end

p = 0*gamma;

dx = abs(X(1,2,1)-X(1,1,1))
de = abs(eta_i(2)-eta_i(1))

a = 0.5*dcx/de;
b = 0.5*dce/dx;

a(isnan(a))=0;
b(isnan(b))=0;
gamma(isnan(gamma))=0;

a(isinf(a))=0;
b(isinf(b))=0;
gamma(isinf(gamma))=0;

ncount=[];
err = [];
p = 0*a;

for k = 1:nt
%     if k>1
%         p(:,:,k) = psi(:,:,k-1);
%     end
    [psi(:,:,k),ncount(k),err(k)] = forward_diff_sor_3(medfilt2(squeeze(a(:,:,k))),medfilt2(squeeze(b(:,:,k))),p(:,:,k),medfilt2(squeeze(gamma(:,:,k))),ne,nx);
end

end

