function [Ec,Ek,dsEc,dsEk,dpicdk,dpikdk,P_c,D_c,D_k,lam_c,lam_k,PP] = model_transport(k,n_x,n_k,eta_x,eta_k,Re,Sc)

[Ek,dsEk,~,lam_k] = model_spectrum(k,n_k,eta_k);

[Ec,dsEc,~,lam_c] = model_spectrum(k,n_x,eta_x);

dk = abs(k(2)-k(1));

a1 = 0.001;
b1 = pi;
a2 = -0.03;
b2 = pi/2;

phik = a1*k+b1;
phic = a2*k+b2;

Ekcphi = (conv(sqrt(Ek(:).*exp(1j*phik(:))),sqrt(Ec(:).*exp(1j*phic(:))),'same')*dk).^2;
Ekkphi = (conv(sqrt(Ek(:).*exp(1j*phik(:))),sqrt(Ek(:).*exp(1j*phik(:))),'same')*dk).^2;

[phikc,Ekc] = cart2pol(real(Ekcphi(:)),imag(Ekcphi(:)));
[phikk,Ekk] = cart2pol(real(Ekkphi(:)),imag(Ekkphi(:)));

dpicdk = k(:).*sqrt(Ekc(:).*Ec(:)).*sin(0.5*(phikc(:)-phic(:)));

pi_c = cumtrapz(dpicdk)*dk;

dpikdk = k(:).*sqrt(Ekk(:).*Ek(:)).*sin(0.5*(phikk(:)-phik(:)));

pi_k = cumtrapz(dpikdk)*dk;

P_c = pi*sqrt(Ek(:).*Ec(:)).*cos(0.5*(phik(:)-phic(:)));

D_c = -1*(k(:).^2).*Ec(:)/(Re*Sc);
D_k = -1*(k(:).^2).*Ek(:)/(Re);

PP = -1*real(sqrt(Ek(:).*exp(1j*phik(:))).*(conv(sqrt((k(:).^2).*Ek(:).*exp(1j*phik(:))),sqrt((k(:).^2).*Ek(:).*exp(1j*phik(:))),'same')*dk)./k(:));

% [F1a,F2a] = F1_F2(n_x,eta_x);
% 
% dpidk = k(:).*sqrt(Ex(:).*Ek(:)) - k(:).*sqrt(Ex(:)).*dsEk(:).*F1a + k(:).*sqrt(Ex(:)).*d2sEk(:).*F2a;
% 
% pi_c = cumtrapz(dpidk);
% 
% [F1a,F2a] = F1_F2(n_k,eta_k);
% 
% dpidk = k(:).*sqrt(Ek(:).*Ek(:)) - k(:).*sqrt(Ek(:)).*dsEk(:).*F1a + k(:).*sqrt(Ex(:)).*d2sEk(:).*F2a;
% 
% pi_k = cumtrapz(dpidk);

end

