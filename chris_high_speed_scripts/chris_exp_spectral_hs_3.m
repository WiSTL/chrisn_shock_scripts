function [Es,D_hat,T_hat,P,Em,Dm,Pm,K] = chris_exp_spectral_hs_3(C,Cm,R,E,T,X,V0,C_h,h0,hd0,Re_hSc,h,l,dpxdx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% version3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nz,nx,nt] = size(C);

for i = 1:nt
    for j = 1:nz
        [Es(j,:,i),K(j,:,i)] = pspectrum(C(j,:,i),(1/0.25)*[1:nx]/dpxdx,'leakage',0.85);
        %[Es(j,:,i),K(j,:,i)] = power_spectra_1D(C(j,:,i),1/dpxdx,1,1);
    end
end

Es = imresize3(Es,[nz,nx,nt]);
K = imresize3(K,[nz,nx,nt]);

Es = medfilt3(Es);

[dee,~,det] = chris_gradient_hs(Es,X,E,T);
[d2ee2,~,~] = chris_gradient_hs(dee,X,E,T);

for i = 1:nt
    detau(:,:,i) = det(:,:,i)*h0/hd0 + dee(:,:,i).*(V0(i)*h0/(h(i)*hd0) + E(:,:,i)*C_h(i));
end

detau = medfilt3(detau);


rhom = 1+(R-1)*Cm;

for k = 1:nt
    for i = 1:nx
        D_hat(:,i,k) = (1/Re_hSc(k))*(1./rhom(:,k)).*(-1*(K(:,i,k).^2).*Es(:,i,k)*(h(k)/l)^2 + d2ee2(:,i,k));
        T_hat(:,i,k) = D_hat(:,i,k) + E(:,i,k).*d2ee2(:,i,k) - (1/C_h(k))*detau(:,i,k);
    end
end

for k = 1:nt
    for j = 1:nz
        P(j,:,k) = cumtrapz(K(j,:,k),-1*T_hat(j,:,k));
    end
end

for i = 1:nt
    for j = 1:nx
        itop = find(E(:,j,i)>-1);
        range_cp = find(E(itop,j,i)<1);

        if ~isempty(range_cp)
            Em(j,i) = -1*trapz(E(range_cp,j,i),Es(range_cp,j,i));
            
            Dm(j,i) = -1*trapz(E(range_cp,j,i),D_hat(range_cp,j,i));
            
            Pm(j,i) = -1*trapz(E(range_cp,j,i),P(range_cp,j,i));
            
        end
    end
end


end

