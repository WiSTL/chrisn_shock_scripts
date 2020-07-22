function [tuu,tuv,tvv,K] = chris_mean_profile_stresses(u,v,r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - taking the x mean of NS Tau_u, Tau_v 
%   and Tau_r appear as forcing of the 
%   mean rho*u, rho*v and rho variables
%   respectivly.
% - applies for 2D compressible fluid
%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm = mean(r');
ru = mean((r.*u)');
rv = mean((r.*v)');
ruu = mean((r.*u.*u)');
ruv = mean((r.*u.*v)');
rvv = mean((r.*v.*v)');

uf = ru./rm;
vf = rv./rm;
uuf = ruu./rm;
uvf = ruv./rm;
vvf = rvv./rm;

tuu = rm.*(uuf - uf.*uf);
tuv = rm.*(uvf - uf.*vf);
tvv = rm.*(vvf - vf.*vf);

Ku = 0.5*tuu;
Kv = 0.5*tvv;

K = (Ku+Kv)./rm;


end

