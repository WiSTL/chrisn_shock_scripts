function [E,E1,F,F1,C_mikaelian,h_dot,eta,K] = chris_exp_growth_rate_pred_hs(C,A,dV0,dpxdx)

[nz,nx] = size(C);

[~,eta,~,delta,y0] = chris_C(C);

h = delta/dpxdx

Cm = mean(C,2);

cp = [];

for i = 1:nx
    cp(:,i) = C(:,i) - Cm;
end

x = [1:nx]/dpxdx;

itop = find(eta>-0.5);
range_cp = find(eta(itop)<0.5);

for j = 1:nz
   [E(j,:),K(j,:)] = pspectrum(cp(j,:),x,'leakage',0.85); 
   [Et(j,:),K(j,:)] = pspectrum(cp(j,:).^2,x,'leakage',0.85); 
end

E1 = -1*trapz(eta,E);
E2 = -1*trapz(eta,Et);

k_h = K(1,:)*h;
G = 1 - exp(-0.7*k_h);
G1 = 1 - exp(-0.49*k_h);

F = G.*E1;
F1 = G1.*E2;

F_i = trapz(k_h,F);
F1_i = trapz(k_h,F1);

C_mikaelian = ((2.85*abs(A)^1.09)./(log((1+A)/(1-A)))).*F_i + ((10.22*abs(A)^1.09)./(log((1+A)/(1-A)))).*F1_i;

h_dot = C_mikaelian.*dV0.*A;

end

