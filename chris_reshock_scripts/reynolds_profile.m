function [t11,t12,t22,up,vp,K] = reynolds_profile(u,v)

up = u - mean2(u);
vp = v - mean2(v);

t11 = mean((up.^2)');
t12 = mean((up.*vp)');
t22 = mean((vp.^2)');

Ku = 0.5*t11;
Kv = 0.5*t22;

K = Ku+Kv;


end

