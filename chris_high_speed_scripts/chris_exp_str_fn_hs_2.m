function [Sx,Sxn,rx,zeta,alpha,zetan,sxnn] = chris_exp_str_fn_hs_2(C,dpxdx,eta,tau,tau_i)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% version4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta_i = -1:0.01:1;
tau_i = tau_i;

[C_i] = chris_interp_hs(C,eta,tau,eta_i,tau_i);

C_i(isnan(C_i))=0;

[ne,nx,nt] = size(C_i);

ft = fittype('a+b*x');

for i = 1:nt
    [Sx(:,:,i),rx(:,:,i)] = chris_no_structure_fn(C_i(:,:,i),2,dpxdx);%,xzn(:,:,i),xrn(:,:,i)
    [dSz(:,:,i),dSr(:,:,i)] = chris_gradient(Sx(:,:,i),1/dpxdx,0.01);
end

Sx = medfilt3(Sx);

k = log(squeeze(rx(1,2:5,1)))';%exp:15:145 sim:55:64

for i = 1:nt
    for j = 1:ne
        a = log(squeeze(Sx(j,2:5,i)))'; %15:145 55:64
        a(isnan(a)) = 0;
        a(isinf(a)) = 0;
        f = fit(k,a, ft,'StartPoint',[1,3]);
        zeta(j,i) = f.b;
        alpha(j,i) = exp(f.a);
    end
end

for n=2:10
    for i = 1:nt
        [Sxn(:,:,i),~] = chris_no_structure_fn(C_i(:,:,i),n,dpxdx);
        a = log(squeeze(nanmean(Sxn(:,2:5,i),1)))';
        a(isnan(a)) = 0;
        a(isinf(a)) = 0;
        f = fit(k,a, ft,'StartPoint',[1,3]);
        zetan(n,i) = f.b;
    end
    sxnn(:,n) = squeeze(nanmean(Sxn(:,:,nt),1));
end

end

