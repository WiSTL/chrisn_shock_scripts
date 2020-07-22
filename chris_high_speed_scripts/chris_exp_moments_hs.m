function [Cm,cp2,cp3,dwc,DC] = chris_exp_moments_hs(C,R,Re_hSc,h,l,dpxdx,eta,tau,tau_i,h0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% version4
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta_i = -2:0.02:2;
tau_i = tau_i;

[C_i] = chris_interp_hs(C,eta,tau,eta_i,tau_i);

C_i(isnan(C_i))=0;

h = interp1(tau,smooth(h),tau_i);
Re_hSc = interp1(tau,smooth(Re_hSc),tau_i);

hh = log(h/h0);

[ne,nx,nt] = size(C_i);

x_h = (1/l)*[1:nx]/dpxdx;

%[taum,etam] = meshgrid(tau_i,eta_i);
[hhm,etam] = meshgrid(hh,eta_i);
[E,~,X,T] = meshgrid_hs(etam,x_h,tau_i,h);

[dce,~,~] = chris_gradient_hs(C_i,X,E,T);
[d2ce2,~,~] = chris_gradient_hs(dce,X,E,T);

Cm = squeeze(mean(C_i,2));
dcem = squeeze(mean(dce,2));
d2ce2m = squeeze(mean(d2ce2,2));

[dcmt,~] = chris_gradient(Cm,0.02,mean(diff(tau_i)));

rhom = 1+(R-1)*Cm;

for i = 1:nt
    dwc(:,i) = d2ce2m(:,i)./(rhom(:,i)*Re_hSc(i)) + eta_i'.*dcem(:,i);
    DC(:,i) = dcmt(:,i);
    for j = 1:nx
       C_p(:,j,i) = C_i(:,j,i) - Cm(:,i); 
    end
end

cp2 = squeeze(mean(C_p.^2,2));
cp3 = squeeze(mean(C_p.^3,2));



end

