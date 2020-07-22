function [C_interp] = chris_interp_hs(C,eta,tau,eta_interp,tau_interp)

nyi = length(eta_interp);
[~,nx,nt] = size(C);
    
C_interp1 = [];
C_interp2 = [];
     
for i = 1:nx
    for j = 1:nt
        C_interp1(:,i,j) = interp1(eta(:,j),squeeze(C(:,i,j)),eta_interp);
    end
    for j = 1:nyi
        C_interp2(j,i,:) = interp1(tau,squeeze(C_interp1(j,i,:)),tau_interp); 
    end
end
     
C_interp = C_interp2;
end

