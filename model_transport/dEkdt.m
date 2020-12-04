function [Ec,Ek,dEc,dEk,lam_c,lam_k,dEkl,dEcl,dnk,dnc,dlnrelk,dlnsc] = dEkdt(k,rel,scl,nc,nk,Re,Sc)

eta_x = 1/sqrt(rel*scl);
eta_k = 1/sqrt(rel);

dk = abs(k(2)-k(1));

[Ec,Ek,dsEc,dsEk,dpicdk,dpikdk,P_c,D_c,D_k,lam_c,lam_k,PP] = model_transport(k,nc,nk,eta_x,eta_k,Re,Sc);

dEk = -1*dpikdk + D_k - PP;
dEc = -1*dpicdk + D_c + P_c;

klam = 2/lam_k;
klamc = 2/lam_c;

dEk(isinf(dEk))=0;
dEc(isinf(dEc))=0;

dEk(isnan(dEk))=0;
dEc(isnan(dEc))=0;

dnkdt = (k(:)./Ek(:)).*chris_derivative(dEk(:),dk) - (k(:)./(Ek(:).^2)).*dEk(:).*dsEk(:);
dncdt = (k(:)./Ec(:)).*chris_derivative(dEc(:),dk) - (k(:)./(Ec(:).^2)).*dEc(:).*dsEc(:);

ii = find(k<klam);
ii2 = find(k(ii)>klam-0.5);

ii3 = ii2(end);

dEkl = dEk(ii3);
dEcl = dEc(ii3);

dnk = dnkdt(ii3);
dnc = dncdt(ii3);


dlamk = trapz(k,dEk)/(2*lam_k*trapz((k.^2).*Ek))*dk - (1/(2*dk*trapz((k.^2).*Ek)))*trapz((k.^2)*dEk)*dk;
dlamc = trapz(k,dEc)/(2*lam_c*trapz((k.^2).*Ec))*dk - (1/(2*dk*trapz((k.^2).*Ec)))*trapz((k.^2)*dEc)*dk;

dlnrelk = dlamk*(0.018*(rel^0.035)*(nk-17.74))^(-1) - (0.0515/0.018)*(1/(nk-17.74))*dnk;

dlnresc = dlamc*(0.018*((rel.*scl)^0.035)*(nc-17.74))^(-1) - (0.0515/0.018)*(1/(nc-17.74))*dnc;

drel = rel*dlnrelk;
drelsc = rel*scl*dlnresc;

dsc = (drelsc/rel) - (scl/rel)*drel;

dlnsc = dsc/scl;

end

