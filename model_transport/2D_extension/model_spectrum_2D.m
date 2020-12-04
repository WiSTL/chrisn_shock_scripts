function [Eu,Ew,Ec,lambdau,lambdaw,lambdac] = model_spectrum_2D(k,z,n_u,eta_L_u,n_w,eta_L_w,n_c,eta_L_c)

[K,Z] = meshgrid(k,z);


Eu = ((0.2+Z.^2).*exp(-1*Z.^2)).*(K.^(-n_u)).*((K./sqrt(K.^2 + 6.78)).^(2+n_u)).*exp(-5.2*eta_L_u.*K);
Eu(isnan(Eu))=0;

Ew = (exp(-1*Z.^2)).*(K.^(-n_w)).*((K./sqrt(K.^2 + 6.78)).^(2+n_w)).*exp(-5.2*eta_L_w.*K);
Ew(isnan(Ew))=0;

Ew1 = Ew/trapz(z,trapz(k,Ew'+Eu'));

Eu = Eu/trapz(z,trapz(k,Eu'+Ew'));

Ew = Ew1;

Ec = (exp(-1*Z.^2)).*(K.^(-n_c)).*((K./sqrt(K.^2 + 6.78)).^(2+n_c)).*exp(-5.2*eta_L_c.*K);
Ec(isnan(Ec))=0;

Ec = Ec/trapz(z,trapz(k,Ec'));

lambdau = nanmean(sqrt(squeeze(trapz(k,Eu'))./squeeze(trapz(k,(K'.^2).*Eu'))));
lambdaw = nanmean(sqrt(trapz(k,Ew')./trapz(k,(K'.^2).*Ew')));
lambdac = nanmean(sqrt(trapz(k,Ec')./trapz(k,(K'.^2).*Ec')));


end