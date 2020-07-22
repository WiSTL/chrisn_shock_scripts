function [w,u,a,c,b,gamma] = chris_velocity_hs_3(C,x,eta,tau,h,L)

eta_i = -2:0.05:2;

[C_i] = chris_interp_hs(C,eta,tau,eta_i,tau);

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
        gamma(j,:,k) = (eta_i(j).*dce(j,:,k) - dct(j,:,k));
    end
end

p = 0*gamma;
u = 0*gamma;
w = 0*gamma;

dx = abs(X(1,2,1)-X(1,1,1))
de = abs(eta_i(2)-eta_i(1))

a = dcx;
b = dce;
c = (b - (dx/de)*a);

a(isnan(a))=0;
b(isnan(b))=0;
c(isnan(c))=0;
gamma(isnan(gamma))=0;

a(isinf(a))=0;
b(isinf(b))=0;
c(isinf(c))=0;
gamma(isinf(gamma))=0;

gamma = medfilt3(gamma);
a = medfilt3(a);
b = medfilt3(b);
c = medfilt3(c);

b = b./c;
a = a./c;
gamma = gamma./c;

a(isnan(a))=0;
b(isnan(b))=0;
gamma(isnan(gamma))=0;

a(isinf(a))=0;
b(isinf(b))=0;
gamma(isinf(gamma))=0;

gamma = medfilt3(gamma(end:-1:1,:,:));
a = medfilt3(a(end:-1:1,:,:));
b = medfilt3(b(end:-1:1,:,:));

for n = 1:nt
    for j = 2:ne
        for i = 2:nx
            w(j,i,n) = (gamma(j,i,n) - w(j-1,i,n)*squeeze(b(j,i,n))*(dx/de) - u(j,i-1,n)*squeeze(a(j,i,n)));
            u(j,i,n) = (w(j-1,i,n) - w(j,i,n))*(dx/de)*(L/h(n)) + u(j,i-1,n);
        end
    end
end

% ncount=[];
% err = [];
% 
% for k = 1:nt
% %     if k>1
% %         p(:,:,k) = psi(:,:,k-1);
% %     end
%     [psi(:,:,k),ncount(k),err(k)] = forward_diff_sor_2(a(:,:,k),b(:,:,k),c(:,:,k),p(:,:,k),de*gamma(:,:,k),ne,nx);
% end

end

