function [Ep,Piu,Piw,Pizu,Pizw,PU,PW,dEu,dEw,dnu,dnw,dnc,dlnrel,dlnsc] = model_transport_terms_2D(k,z,n_u,eta_L_u,n_w,eta_L_w,n_c,eta_L_c,R,Re,Sc)

[Eu,Ew,Ec,lam_u,lam_w,lam_c] = model_spectrum_2D(k,z,n_u,eta_L_u,n_w,eta_L_w,n_c,eta_L_c);

[K,Z] = meshgrid(k,z);

dk = abs(k(2)-k(1));
dz = abs(z(2)-z(1));

[dEz,dEk] = chris_gradient(sqrt(Ew),dz,dz);

for i = 1:length(z)
    Ep(i,:) = (1+0.5*(R-1)*erfc(z(i))).*(conv(sqrt(squeeze(Eu(i,:))),sqrt(squeeze(Eu(i,:))),'same')*dk - conv(squeeze(dEz(i,:)),squeeze(dEz(i,:)),'same').*dk./(k.^2)).^2;
    
    Ep(isnan(Ep))=0;
    Ep(isinf(Ep))=0;
    
    Euu(i,:) = (conv(sqrt(squeeze(Eu(i,:))),sqrt(squeeze(Eu(i,:))),'same')*dk).^2;
    Euw(i,:) = (conv(sqrt(squeeze(Eu(i,:))),sqrt(squeeze(Ew(i,:))),'same')*dk).^2;
    Ewu(i,:) = (conv(sqrt(squeeze(Ew(i,:))),sqrt(squeeze(Eu(i,:))),'same')*dk).^2;
    Eww(i,:) = (conv(sqrt(squeeze(Ew(i,:))),sqrt(squeeze(Ew(i,:))),'same')*dk).^2;
    
    Ewc(i,:) = (conv(sqrt(squeeze(Ew(i,:))),sqrt(squeeze(Ec(i,:))),'same')*dk).^2;
    Euc(i,:) = (conv(sqrt(squeeze(Eu(i,:))),sqrt(squeeze(Ec(i,:))),'same')*dk).^2;
    
    Epu(i,:) = (conv(sqrt(squeeze(Ep(i,:))),sqrt(squeeze(Eu(i,:))),'same')*dk).^2;
    Epw(i,:) = (conv(sqrt(squeeze(Ep(i,:))),sqrt(squeeze(Ew(i,:))),'same')*dk).^2;
end

dPiu = K.*sqrt(Euu.*Eu);
dPiw = K.*sqrt(Euw.*Ew);
dPic = K.*sqrt(Euc.*Ec);

Piu = cumtrapz(k,dPiu')';
Piw = cumtrapz(k,dPiw')';
Pic = cumtrapz(k,dPic')';

Piu(isnan(Piu))=0;
Piw(isnan(Piw))=0;
Pic(isnan(Pic))=0;

Pizu = sqrt(Ewu.*Eu);
Pizw = sqrt(Eww.*Ew);
Pizc = sqrt(Ewc.*Ec);

Pizu(isnan(Pizu))=0;
Pizw(isnan(Pizw))=0;
Pizc(isnan(Pizc))=0;

Ep(isnan(Ep))=0;
Ep(isinf(Ep))=0;

PU = K.*sqrt(Ep.*Eu)./(1+0.5*(R-1)*erfc(Z));
PW = sqrt(Ep.*Ew);

PW = chris_gradient(PW,dz,dz)./(1+0.5*(R-1)*erfc(Z));

Tu = dPiu;
Tw = dPiw;
Tc = dPic;

[Tzu,~] = chris_gradient(Pizu,dz,dz);
[Tzw,~] = chris_gradient(Pizu,dz,dz);
[Tzc,~] = chris_gradient(Pizc,dz,dz);

Du = -1*(K.^2).*Eu./((1+0.5*(R-1)*erfc(Z))*Re);
Dw = -1*(K.^2).*Ew./((1+0.5*(R-1)*erfc(Z))*Re);
Dc = -1*(K.^2).*Ec./(Sc*Re);

[dEudz,~] = chris_gradient(Eu,dz,dz);
[dEwdz,~] = chris_gradient(Ew,dz,dz);
[dEcdz,~] = chris_gradient(Ec,dz,dz);

Dzu = chris_gradient(dEudz,dz,dz)./((1+0.5*(R-1)*erfc(Z))*Re);
Dzw = chris_gradient(dEwdz,dz,dz)./((1+0.5*(R-1)*erfc(Z))*Re);
Dzc = chris_gradient(dEcdz,dz,dz)./(Sc*Re);

[dsEu,~] = chris_gradient(sqrt(Eu),dz,dz);
[dsEw,~] = chris_gradient(sqrt(Ew),dz,dz);
[dsEc,~] = chris_gradient(sqrt(Ec),dz,dz);

Xu = 2*(dsEu.^2)./((1+0.5*(R-1)*erfc(Z))*Re);
Xw = 2*(dsEw.^2)./((1+0.5*(R-1)*erfc(Z))*Re);
Xc = 2*(dsEc.^2)./(Sc*Re);

[dCm,~] = chris_gradient(0.5*erfc(Z),dz,dz);
Pc = sqrt(Ec.*Ew).*dCm;

dEu = Tu + Tzu + Du + Dzu - PU - Xu;
dEw = Tw + Tzw + Dw + Dzw - PW - Xw;
dEc = Pc + Tc + Tzc + Dc + Dzc - Xc;

[~,d2Eudk] = chris_gradient(dEu,dk,dk);
[~,dEudk] = chris_gradient(Eu,dk,dk);

[~,d2Ewdk] = chris_gradient(dEw,dk,dk);
[~,dEwdk] = chris_gradient(Ew,dk,dk);

[~,d2Ecdk] = chris_gradient(dEc,dk,dk);
[~,dEcdk] = chris_gradient(Ec,dk,dk);


dnu1 = (K./Eu).*d2Eudk - (K./(Eu.^2)).*dEu.*dEudk;
dnw1 = (K./Ew).*d2Ewdk - (K./(Ew.^2)).*dEw.*dEwdk;
dnc1 = (K./Ec).*d2Ecdk - (K./(Ec.^2)).*dEc.*dEcdk;

klamu = 2/lam_u;
klamw = 2/lam_w;
klamc = 2/lam_c;

ii = find(k<klamu);
ii2 = find(k(ii)>klamu-0.5);

ii3u = ii2(end);

ii = find(k<klamw);
ii2 = find(k(ii)>klamw-0.5);

ii3w = ii2(end);

ii = find(k<klamc);
ii2 = find(k(ii)>klamc-0.5);

ii3c = ii2(end);

for i = 1:length(z)
    dnu(i) = dnu1(i,ii3u);
    dnw(i) = dnw1(i,ii3w);
    dnc(i) = dnc1(i,ii3c);
end

rel = 1./(eta_L_w.^2);
relsc = 1./(eta_L_c.^2);

scl = relsc/rel;

dlamw = trapz(k,dEw')./(2*lam_w*trapz(((K.^2).*Ew)'))*dk - (1./(2*dk*trapz(((K.^2).*Ew)'))).*trapz(((K.^2).*dEw)')*dk;
dlamc = trapz(k,dEc')./(2*lam_c*trapz(((K.^2).*Ec)'))*dk - (1./(2*dk*trapz(((K.^2).*Ec)'))).*trapz(((K.^2).*dEc)')*dk;

dlnrel = dlamw*(0.018*(rel^0.035)*(n_w-17.74))^(-1) - (0.0515/0.018)*(1/(n_w-17.74))*dnw;
dlnrelsc = dlamc*(0.018*(relsc^0.035)*(n_c-17.74))^(-1) - (0.0515/0.018)*(1/(n_c-17.74))*dnc;

drel = rel*dlnrel;
drelsc = rel*scl*dlnrelsc;

dsc = (drelsc/rel) - (scl/rel)*drel;

dlnsc = dsc/scl;

end

