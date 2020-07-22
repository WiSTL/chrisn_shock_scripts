function [Em,D_hat,T_hat,P,kx] = chris_exp_spectral_hs(C,rho_bar,eta,E,T,X,C_h,h0,hd0,Re_hSc,hl,dpxdx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% version1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nz,nx,nt] = size(C);

for i = 1:nt
    for j = 1:nz
        [Es(j,:,i),K(j,:,i)] = pspectrum(C(j,:,i),(1/0.25)*[1:nx]/dpxdx);
    end
end

Es = imresize3(Es,[nz,nx,nt]);
K = imresize3(K,[nz,nx,nt]);

[dee,~,~] = chris_gradient_hs(Es,X,E,T);
[d2ee2,~,~] = chris_gradient_hs(dee,X,E,T);


for i = 1:nt
    itop = find(eta(:,i)>-0.5);
    range_cp = find(eta(itop,i)<0.5);

    if ~isempty(range_cp)
        for j = 1:nx
            Em(j,i) = trapz(eta(range_cp,i),Es(range_cp,j,i));
            Em2(j,i) = trapz(eta(range_cp,i),Es(range_cp,j,i)./rho_bar(range_cp,i));
            Em3(j,i) = trapz(eta(range_cp,i),d2ee2(range_cp,j,i)./rho_bar(range_cp,i));
            Em4(j,i) = trapz(eta(range_cp,i),eta(range_cp,i).*dee(range_cp,j,i));
        end
        kx(:,i) = mean(K(range_cp,:,i),1);
    end
end

[dest,~] = chris_gradient(Em,kx,squeeze(T(1,:,:)));

destau = dest*h0/hd0;

for i = 1:nt
    D_hat(:,i) = (hl(i)/Re_hSc(i))*((kx(:,i).^2).*Em2(:,i));
    T_hat(:,i) = D_hat(:,i) - Em(:,i) - (1/C_h(i))*destau(:,i);
end

for i = 1:nt
    P(:,i) = cumtrapz(kx(:,i),T_hat(:,i));
end

end

