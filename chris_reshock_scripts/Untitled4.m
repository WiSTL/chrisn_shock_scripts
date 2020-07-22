[up2n,upvpn,apupn,apvpn,vp2n,vp3n,upupvpn,upvpdvpn,vpvpdupn,cn,cpn,un,vn,xin,udcn,vdcn,t11n,t12n,t22n,bn,t11Rn,t12Rn,t22Rn,rvn,run,eresn,esgsn,b11Rn,b12Rn,b22Rn,q2n,eta_m,Re] = chris_exp_profiles(1,[-0.8:0.005:0.8]);

eta = eta_m;

vn = mean(vn);
un = mean(un);
cn = mean(cn);

upvpn = mean(upvpn);
up2n = mean(up2n);
vp2n = mean(vp2n);

vp3n = mean(vp3n);

upupvpn = mean(upupvpn);

[dvm] = chris_derivative(smooth(vn,20)',abs(eta(2)-eta(1)));
[dum] = chris_derivative(smooth(un,20)',abs(eta(2)-eta(1)));
[dcm] = chris_derivative(smooth(cn,20)',abs(eta(2)-eta(1)));

[d2vm] = chris_derivative(smooth(dvm,20)',abs(eta(2)-eta(1)));
[d2um] = chris_derivative(smooth(dum,20)',abs(eta(2)-eta(1)));
[d2cm] = chris_derivative(smooth(dcm,20)',abs(eta(2)-eta(1)));

tkk = up2n+vp2n;

dT = chris_derivative(smooth(upupvpn+(1/3)*vp3n,20)',abs(eta(2)-eta(1)));
