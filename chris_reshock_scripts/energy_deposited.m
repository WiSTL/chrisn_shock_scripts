function [vp,dvp,lambda,f1,f2,eta1,eta2] = energy_deposited(kh,R)

eta2 = -10:0.001:0;
eta1 = 0:0.001:10;

eta3 = -10:0.001:10;

f1 = kh*exp(-kh*eta1).*(abs(R-1)*exp(-eta1.^2))./((1+0.5*(R-1).*erfc(eta1)).^1);
f2 = kh*exp(kh*eta2).*(abs(R-1)*exp(-eta2.^2))./((1+0.5*(R-1).*erfc(eta1)).^1);

f3 = -(1./kh).*exp(-kh*eta1).*(abs(R-1).*eta1.*exp(-eta1.^2))./((1+0.5*(R-1).*erfc(eta1)).^2);
f4 = -(1./kh).*exp(kh*eta2).*(abs(R-1).*eta2.*exp(-eta2.^2))./((1+0.5*(R-1).*erfc(eta2)).^2);

f5 = -(1./kh).*exp(-kh*eta1).*(((R-1).^2).*exp(-2*eta1.^2))./((1+0.5*(R-1).*erfc(eta1)).^3);
f6 = -(1./kh).*exp(kh*eta2).*(((R-1).^2).*exp(-2*eta2.^2))./((1+0.5*(R-1).*erfc(eta2)).^3);

f7 = - (1./kh).*(((R-1).^2))./((1+0.5*(R-1)).^3);

g1 = (kh^2)*exp(-kh*eta1).*(abs(R-1)*exp(-eta1.^2))./((1+0.5*(R-1).*erfc(eta1)).^2);
g2 = (kh^2)*exp(kh*eta2).*(abs(R-1)*exp(-eta2.^2))./((1+0.5*(R-1).*erfc(eta1)).^2);



vp = trapz(eta1,f1) + trapz(eta2,f2);% + trapz(eta1,f3) + trapz(eta2,f4) + trapz(eta1,f5) + trapz(eta2,f6) + f7;
dvp = sqrt(trapz(eta1,g1).^2 + trapz(eta2,g2).^2);

lambda = vp./dvp;

end

