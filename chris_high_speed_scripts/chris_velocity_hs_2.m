function [psi,ncount,err,a,c,b,gamma] = chris_velocity_hs_2(C,x,eta,tau,Ch,delta,hl)

eta_i = -2:0.05:2;
tau_i = 0:0.05:5;

[C_i] = chris_interp_hs(C,eta,tau,eta_i,tau_i);

Ch_i = interp1(tau,smooth(Ch),tau_i);
hl_i = interp1(tau,smooth(hl),tau_i);
delta_i = interp1(tau,smooth(delta),tau_i);

[ne,nx,nt] = size(C_i);

[~,etam] = meshgrid(tau_i,eta_i);

[E,Z,X,T] = meshgrid_hs(etam,x,tau_i,delta_i);

[dce,dcx,dct] = chris_gradient_hs(C_i,X,E,T);

dce = medfilt3(dce);
dcx = medfilt3(dcx);
dct = medfilt3(dct);



for k = 1:nt
    for j = 1:ne
        gamma(j,:,k) = (eta_i(j).*dce(j,:,k) - (1/Ch_i(k)).*dct(j,:,k))./hl_i(k);
    end
end

p = 0*gamma;

dx = abs(X(1,2,1)-X(1,1,1))
de = abs(eta_i(2)-eta_i(1))

a = dcx/de;
c = dce/dx;
b = a+c;

a(isnan(a))=0;
b(isnan(b))=0;
c(isnan(c))=0;
gamma(isnan(gamma))=0;

a(isinf(a))=0;
b(isinf(b))=0;
c(isinf(c))=0;
gamma(isinf(gamma))=0;

ncount=[];
err = [];

for k = 1:nt
%     if k>1
%         p(:,:,k) = psi(:,:,k-1);
%     end
    [psi(:,:,k),ncount(k),err(k)] = forward_diff_sor(a(:,:,k),b(:,:,k),c(:,:,k),p(:,:,k),gamma(:,:,k),ne,nx);
end

end

