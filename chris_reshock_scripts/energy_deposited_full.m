function [vp,eta1,eta2] = energy_deposited_full(k,h,R)

eta2 = -10:0.001:0;
eta1 = 0:0.001:10;

rhoh1 = abs(R-1)*exp(-((eta1 - 0.32)/0.25).^2)*k^-0.64;
rhoh2 = abs(R-1)*exp(-((eta2 - 0.32)/0.25).^2)*k^-0.64;

rhom1 = 1+0.5*(R-1).*erfc(eta1);
rhom2 = 1+0.5*(R-1).*erfc(eta2);

f1 = k*h*exp(-k*h*eta1).*rhoh1./(rhom1).^2;
f2 = k*h*exp(-k*h*eta2).*rhoh2./(rhom2).^2;

vp = trapz(eta1,f1) + trapz(eta2,f2);

end

