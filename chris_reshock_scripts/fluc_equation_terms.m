function [P,T,D,e,dd,up2,up3,vp2,vp3,cp2,cp3,upvp] = fluc_equation_terms(u,v,C,rho,deltadot,delta,nu,eta)

[ny,nx] = size(u);

u = u/deltadot;
v = v/deltadot;

rho0 = 7.622;

rho = rho/rho0;
alpha = 1./rho;

Re = delta*deltadot/mean2(nu);

um = mean(u');
vm = mean(v');
Cm = mean(C');
am = mean(alpha');
rm = mean(rho');

for j = 1:nx
u_p(:,j) = u(:,j) - um';
end

for j = 1:nx
v_p(:,j) = v(:,j) - vm';
end

for j = 1:nx
C_p(:,j) = C(:,j) - Cm';
end

for j = 1:nx
a_p(:,j) = alpha(:,j) - am';
end

for j = 1:nx
r_p(:,j) = rho(:,j) - rm';
end

cp2 = mean((medfilt2(C_p.^2))');
cp3 = mean((medfilt2(C_p.^3))');
cpvp = mean((medfilt2(C_p.*v_p))');
apvp = mean((medfilt2(a_p.*v_p))');
apup = mean((medfilt2(a_p.*u_p))');
up2 = mean((medfilt2(u_p.^2))');
up3 = mean((medfilt2(u_p.^3))');
vp2 = mean((medfilt2(v_p.^2))');
upvp = mean((medfilt2(u_p.*v_p))');
vp3 = mean((medfilt2(v_p.^3))');
rpup = mean((medfilt2(r_p.*u_p))');
rpvp = mean((medfilt2(r_p.*v_p))');
rpup2 = mean((medfilt2(r_p.*u_p.^2))');
rpvp2 = mean((medfilt2(r_p.*v_p.^2))');
rpup3 = mean((medfilt2(r_p.*u_p.^3))');
rpvp3 = mean((medfilt2(r_p.*v_p.^3))');


upupvp = mean((medfilt2(u_p.*u_p.*v_p))');
upvpvp = mean((medfilt2(u_p.*v_p.*v_p))');
cpvpvp = mean((medfilt2(C_p.*v_p.*v_p))');

[drpvp] = chris_derivative(smooth(rpvp)',abs(eta(2)-eta(1)));

[dvp2] = chris_derivative(smooth(vp2)',abs(eta(2)-eta(1)));

[dvp,~] = chris_gradient(medfilt2(v_p),abs(eta(2)-eta(1)),abs(eta(2)-eta(1)));
[dup,~] = chris_gradient(medfilt2(u_p),abs(eta(2)-eta(1)),abs(eta(2)-eta(1)));
[dcp,~] = chris_gradient(medfilt2(C_p),abs(eta(2)-eta(1)),abs(eta(2)-eta(1)));

[dvm] = chris_derivative(smooth(vm)',abs(eta(2)-eta(1)));
[dum] = chris_derivative(smooth(um)',abs(eta(2)-eta(1)));
[dcm] = chris_derivative(smooth(Cm,100)',abs(eta(2)-eta(1)));
[dam] = chris_derivative(smooth(am,100)',abs(eta(2)-eta(1)));

[dvm2] = chris_derivative(smooth(vm.^2)',abs(eta(2)-eta(1)));
[d2vm2] = chris_derivative(smooth(dvm2)',abs(eta(2)-eta(1)));

[d2vm] = chris_derivative(smooth(dvm)',abs(eta(2)-eta(1)));
[d2um] = chris_derivative(smooth(dum)',abs(eta(2)-eta(1)));
[d2cm] = chris_derivative(smooth(dcm,100)',abs(eta(2)-eta(1)));

duuv = chris_derivative(upupvp,abs(eta(2)-eta(1)));
duvv = chris_derivative(upvpvp,abs(eta(2)-eta(1)));
dcvv = chris_derivative(cpvpvp,abs(eta(2)-eta(1)));
dvvv = chris_derivative(vp3,abs(eta(2)-eta(1)));
dcc = chris_derivative(0.5*cp2,abs(eta(2)-eta(1)));
d2uu = chris_derivative(chris_derivative(up2,abs(eta(2)-eta(1))),abs(eta(2)-eta(1)));
d2uv = chris_derivative(chris_derivative(mean((u_p.*v_p)'),abs(eta(2)-eta(1))),abs(eta(2)-eta(1)));
d2cv = chris_derivative(chris_derivative(mean((C_p.*v_p)'),abs(eta(2)-eta(1))),abs(eta(2)-eta(1)));
d2vv = chris_derivative(chris_derivative(vp2,abs(eta(2)-eta(1))),abs(eta(2)-eta(1)));

[d2up,~] = chris_gradient(medfilt2(dup),abs(eta(2)-eta(1)),abs(eta(2)-eta(1)));
[d2vp,~] = chris_gradient(medfilt2(dvp),abs(eta(2)-eta(1)),abs(eta(2)-eta(1)));

upupdvp = mean((u_p.*u_p.*dvp)');
upvpdvp = mean((u_p.*v_p.*dvp)');
vpvpdup = mean((v_p.*v_p.*dup)');
cpvpdcp = mean((C_p.*v_p.*dcp)');
cpvpdvp = mean((C_p.*v_p.*dvp)');
vpdcp = mean((v_p.*dcp)');

apd2up = mean((a_p.*d2up)');
apd2vp = mean((a_p.*d2vp)');


P(:,1) = -0.5*vm.*dvp2';
T(:,1) = vm.^2;
D(:,1) = 0.5*d2vm2'./(Re.*rm);
e(:,1) = (dvm'.^2)./(Re.*rm);

P(:,2) = -0.5*vp2.*dvm';
T(:,2) = vp2;
D(:,2) = 0.5*d2vm2'./(Re.*rm);
e(:,2) = mean((dvp.^2)')./(Re.*rm);

P(:,3) = -dvm'./rm;
T(:,3) = 1./(rm*(0.4));
D(:,3) = drpvp'./rm.^2;
e(:,3) = mean((dvp.^2)')./(Re.*rm);

% P(:,1) = -1*upvp.*dum;
% T(:,1) = -1*0.5*duuv + 0.5*upupdvp;
% D(:,1) = (d2uu);
% e(:,1) = mean((dup.^2)');
% 
% P(:,2) = -1*(up2+vp2).*dum - upvp.*dvm;
% T(:,2) = -1*0.5*duuv;
% D(:,2) = (d2uv);
% e(:,2) = 2*mean((dup.*dvp)');
% 
% P(:,3) = -1*vp2.*dvm;
% T(:,3) = -1*0.5*dvvv;
% D(:,3) = (d2vv);
% e(:,3) = 2*mean((dvp.^2)');

P(:,4) = -1*cpvp.*dcm';
T(:,4) = -1*cpvpdcp;
D(:,4) = vpdcp;
e(:,4) = vm.*dcm';

P(:,5) = -1*vp2.*dcm' - cpvp.*dvm';
T(:,5) = -1*dcvv' + cpvpdvp;
D(:,5) = d2cv';
e(:,5) = 2*mean((dcp.*dvp)');

tkk = up2+vp2;

dtkk = chris_derivative(smooth(tkk,100)',abs(eta(2)-eta(1)));


dd = -(vm.*(dcc)*delta) - delta*(cpvp.*dcm+cpvpdcp);
term = mean((r_p.^2)')./rm.^2;
eps = dam.*dtkk;

end

