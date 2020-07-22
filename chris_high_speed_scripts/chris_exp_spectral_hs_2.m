function [Em,D_hat,T_hat,P,kx] = chris_exp_spectral_hs_2(C,R,eta,E,T,X,C_h,h0,hd0,Re_hSc,hl,dpxdx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% version2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nz,nx,nt] = size(C);
[dce,dcx,~] = chris_gradient_hs(C,X/0.25,E,T);
for i = 1:nt
    a(:,:,i) = (1./Re_hSc(i)).*((R-1)./(1+(R-1.*C(:,:,i)))).*((dcx(:,:,i).^2).*(hl(i).^2)+(dce(:,:,i).^2));
end

for i = 1:nt
    for j = 1:nz
        %[Es(j,:,i),K(j,:,i)] = pspectrum(C(j,:,i),(1/0.25)*[1:nx]/dpxdx,'leakage',0.85);
        [Es(j,:,i),K(j,:,i)] = power_spectra_1D(C(j,:,i),1/dpxdx,1,1);
        a_hat =  fft_corrected(a(j,:,i),1);
        c_hat_s =  conj(fft_corrected(C(j,:,i),1));
        
        Ed(j,:,i) = c_hat_s.*a_hat;
    end
end

Es = imresize3(Es,[nz,nx,nt]);
K = imresize3(K,[nz,nx,nt]);

[dee,~,~] = chris_gradient_hs(Es,X,E,T);

for i = 1:nt
    itop = find(eta(:,i)>-1);
    range_cp = find(eta(itop,i)<1);

    if ~isempty(range_cp)
        for j = 1:nx
            Em(j,i) = trapz(eta(range_cp,i),Es(range_cp,j,i));
            Em2(j,i) = trapz(eta(range_cp,i),eta(range_cp,i).*dee(range_cp,j,i));
            Em3(j,i) = trapz(eta(range_cp,i),Ed(range_cp,j,i));
        end
        kx(:,i) = mean(K(range_cp,:,i),1);
    end
end

[dest,~] = chris_gradient(Em,kx,squeeze(T(1,:,:)));

destau = dest*h0/hd0;

for i = 1:nt
    D_hat(:,i) = -1*(hl(i)/Re_hSc(i))*((kx(:,i).^2).*Em(:,i)) ;%- abs(Em3(:,i));
    T_hat(:,i) = D_hat(:,i) - Em2(:,i) - (1/C_h(i))*destau(:,i);
end

for i = 1:nt
    P(:,i) = cumtrapz(kx(:,i),-1*T_hat(:,i))-1*T_hat(1,i);
end

end

