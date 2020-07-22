function [S, Smax, R] = C_shannon_entropy(C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Shannon Information Entropy
%%% Taken from 
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx,ny] = size(C);

C2 = 1-C;

countzero = length(C(C==0));
countone = length(C(C==1));

% eta = (C.*C2)/0.25;
% eta = sum(sum(subplus(eta)))/(nx*ny-(countone+countzero));
[Cm1, Cm2] = mixing_layer_mask(C,C2);
%eta = var(subplus(Cm2(:)))/mean(mean(subplus(Cm2)));%mixing_measure(subplus(C.*C2));

% rho0 = 1.784;
% rho1 = 0.1786;
% rho = rho0+(rho1-rho0)*C;
% [drdx,drdy] = chris_gradient(rho,1,1);
% eta = sum(sum(abs(medfilt2(drdy,[5 5]))./rho.^2));

C(C<eps)=eps;
C2(C2<eps)=eps;

S=-sum(sum(C.*log(C) + C2.*log(C2) ))/(nx*ny-(countone+countzero))/log(2);

Cm1 = sum(sum(C))/(nx*ny-(countone+countzero));
Cm2 = sum(sum(1-C))/(nx*ny-(countone+countzero));

Smax = -(Cm1*log(Cm1)+Cm2*log(Cm2))/log(2);

% Smax = sum(sum((rho.^2).*C.*(1-C)))


R = sqrt(sum(sum(dot(gradient(C.^0.5),gradient(C.^0.5)))));



end

