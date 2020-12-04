function [E,dsE,d2sE,lambda,E0] = model_spectrum(k,n,eta_L)

E = (k.^(-n)).*((k./sqrt(k.^2 + 6.78)).^(2+n)).*exp(-5.2*eta_L.*k);

E(isnan(E))=0;

E0 = trapz(k,E);

E = E/trapz(k,E);

lambda = sqrt(trapz(k,E)./trapz(k,(k.^2).*E));

dsE = chris_derivative(E,abs(k(2)-k(1)));
d2sE = chris_derivative(dsE,abs(k(2)-k(1)));

end