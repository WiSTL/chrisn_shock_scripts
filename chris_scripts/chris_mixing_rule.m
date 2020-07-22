function [mu,mueff,D,Deff,Sc_eff] = chris_mixing_rule(P,T,xi,xi_acetone)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   1 - Argon
%   2 - Helium
%   3 - Acetone
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M(1) = 39.948;
M(2) = 4.003;
M(3) = 58.0805;

R = 8314;

acetone = xi_acetone;
helium = 1 - acetone;

yac = acetone*M(3)/(acetone*M(3) + helium*M(2));
yhelium = helium*M(2)/(acetone*M(3) + helium*M(2));

mw_mix = 1/((yac/M(3))+(yhelium/M(2)));

e(1) = 119.8;
e(2) = 10.9;
e(3) = 458;

sig(1) = 3.405;
sig(2) = 2.64;
sig(3) = 4.599;

g(1) = 1.667;
g(2) = 1.667;
g(3) = 1.1;

Am = 1.16145; Bm = -0.14874; Cm = 0.52487; Dm = -0.7732; Em = 2.16178;
Fm = -2.43787;

for i = 1:3
    Ts = T/e(i);
    Sig_mu(i) = Am*Ts^Bm + Cm*exp(Dm*Ts) + Em*exp(Fm*Ts);
end

Ad = 1.06036; Bd = -0.1561; Cd = 0.193; Dd = -0.47635; Ed = 1.03587; Fd = -1.52996; 
Gd = 1.76474; Hd = -3.89411;

for i = 1:3
    for j = 1:3
        Ts = T/sqrt(e(i)*e(j));
        Sig_D(i,j) = Ad*Ts^Bd + Cd*exp(Dd*Ts) + Ed*exp(Fd*Ts)+ Gd*exp(Hd*Ts);
    end
end

for i = 1:3
    for j = 1:3
        Mij(i,j) = 2/((1/M(i))+(1/M(j)));
    end
end

for i = 1:3
    mu(i) = (2.6693*10^(-6))*sqrt(M(i)*T)/(Sig_mu(i)*sig(i)^2);
end

for i = 1:3
    for j = 1:3
        sigij = 0.5*(sig(i)+sig(j));
        D(i,j) =   (0.0266/Sig_D(i,j))*(T^(3/2))/(P*sqrt(Mij(i,j))*sigij^2); 
    end
end


Deff = xi./((helium*xi./D(2,1))+(acetone*xi./D(3,1)));

mu_mix = ((mu(3)*yac./M(3)) + (mu(2)*yhelium./M(2)))./((yac/M(3))+(yhelium/M(2)));

y_mix = xi*mw_mix/(xi*mw_mix + (1-xi)*M(1));
y1 = (1-xi)*M(1)/(xi*mw_mix + (1-xi)*M(1));

mueff = ((mu_mix*y_mix./mw_mix)+mu(1)*y1./M(1))./((y_mix./mw_mix)+y1./M(1));

%%%%%%%%%%%%%%%%%%%%%%%%%

yach = (xi)*mw_mix/(xi*mw_mix + (1-xi)*M(1));
yar = (1-xi)*M(1)/(xi*mw_mix + (1-xi)*M(1));

mw_mix_2 = 1/((yach/mw_mix)+(yar/M(1)));

rho = (P/((R/mw_mix_2)*T));

Sc_eff = mueff/(rho*Deff);

end

