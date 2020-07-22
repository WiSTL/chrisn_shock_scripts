function [Es,X_hat,D_hat,T_hat,P,Xm,Em,Dm,Tm,Pm,Xt,Tbt,Ttt,Pt,Dt,K,zeta,alpha,C_i,L] = chris_exp_spectral_hs_4(C,R,Re_hSc,h,l,dpxdx,eta,tau,tau_i,h0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% version4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta_i = -1:0.005:1;
tau_i = tau_i;

[C_i] = chris_interp_hs(C,eta,tau,eta_i,tau_i);

C_i(isnan(C_i))=0;

h = interp1(tau,smooth(h),tau_i);
Re_hSc = interp1(tau,smooth(Re_hSc),tau_i);

lh = log(h/h0);

[ne,nx,nt] = size(C_i);

x_h = (1/l)*[1:nx]/dpxdx;

ft = fittype('a+b*x');

leak = 1;

for i = 1:nt
    [daz,dax] = chris_gradient(squeeze(C_i(:,:,i)),0.005,0.005);
    for j = 1:ne
        [Es(j,:,i),K(j,:,i)] = pspectrum(C_i(j,:,i),x_h,'leakage',leak);
        [Ex(j,:,i)] = pspectrum(dax(j,:),x_h,'leakage',leak);
        [Ez(j,:,i)] = pspectrum(daz(j,:),x_h,'leakage',leak);
        L(j,i) = trapz(squeeze(K(j,2:end,i)),squeeze(Es(j,2:end,i))./squeeze(K(j,2:end,i)))./trapz(squeeze(K(j,2:end,i)),squeeze(Es(j,2:end,i)));
        Xi(j,i) = trapz(squeeze(K(j,2:end,i)),squeeze(K(j,2:end,i).^2).*squeeze(Es(j,2:end,i)));
    end
end

% Es = imresize3(Es,[ne,nx,nt]);
% Ex = imresize3(Ex,[ne,nx,nt]);
% Ez = imresize3(Ez,[ne,nx,nt]);
% K = imresize3(K,[ne,nx,nt]);

[~,nk,~] = size(K);

Es = medfilt3(Es);

[~,etam] = meshgrid(tau_i,eta_i);
[E,~,X,T] = meshgrid_hs(etam,squeeze(K(1,:,1)),lh,h);%tau_i

[dee,~,det] = chris_gradient_hs(Es,X,E,T);
[d2ee2,~,~] = chris_gradient_hs(dee,X,E,T);

detau = medfilt3(det);

for k = 1:nt
    for i = 1:nk
        X_hat(:,i,k) = (1/Re_hSc(k)).*Ez(:,i,k);
        D_hat(:,i,k) = (1/Re_hSc(k)).*(-1*(K(:,i,k).^2).*Es(:,i,k)*(h(k)/l)^2 + d2ee2(:,i,k));
        T_hat(:,i,k) = D_hat(:,i,k) - X_hat(:,i,k) + 0.5*E(:,i,k).*dee(:,i,k) - detau(:,i,k);
    end
end

% for i = 1:nt
%     for j = 1:ne
%         Xi(j,i) = trapz(squeeze(K(j,2:end,i)),squeeze(D_hat(j,2:end,i)));
%     end
% end

for k = 1:nt
    for j = 1:ne
        P(j,:,k) = cumtrapz(K(j,:,k),-1*T_hat(j,:,k));
    end
end

for i = 1:nt
    for j = 1:nk
        itop = find(E(:,j,i)>-1);
        range_cp = find(E(itop,j,i)<1);

        if ~isempty(range_cp)
            Em(j,i) = trapz(E(range_cp,j,i),Es(range_cp,j,i));
            
            Dm(j,i) = -1*trapz(E(range_cp,j,i),D_hat(range_cp,j,i));
            
            Tm(j,i) = -1*trapz(E(range_cp,j,i),T_hat(range_cp,j,i));
            
            Pm(j,i) = -1*trapz(E(range_cp,j,i),P(range_cp,j,i));
            
            Xm(j,i) = -1*trapz(E(range_cp,j,i),X_hat(range_cp,j,i));
            
        end
    end
end

Pt=[];
Dt=[];

for i = 1:nk
    for j = 1:ne
        Pt(j,i) = trapz(lh,squeeze(P(j,i,:)));%trapz(tau_i,C_h'.*squeeze(P(j,i,:)));
        Dt(j,i) = trapz(lh,squeeze(D_hat(j,i,:)));%trapz(tau_i,C_h'.*squeeze(D_hat(j,i,:)));
        Tbt(j,i) = trapz(lh,0.5*squeeze(E(j,i,:).*dee(j,i,:)));
        Ttt(j,i) = trapz(lh,squeeze(T_hat(j,i,:)));
        Xt(j,i) = trapz(lh,squeeze(X_hat(j,i,:)));
    end
end

itop = find(eta_i>-1);
range_cp = find(eta_i(itop)<1);
range_cp = itop(1):itop(range_cp(end));

dk = 1500:4000;%dk = 1500:4000; %exp:15:145 sim:55:64

k = log(squeeze(K(1,dk,1)))';

for i = 1:nt
    for j = 1:ne
        a = log(squeeze(Es(j,dk,i)))'; %15:145 55:64
        a(isnan(a)) = 0;
        a(isinf(a)) = 0;
        f = fit(k,a, ft,'StartPoint',[1,3]);
        zeta(j,i) = f.b;
        alpha(j,i) = exp(f.a);
    end
end

end

