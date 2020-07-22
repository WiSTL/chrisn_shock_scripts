function [G] = chris_zhou_mixed_mass_fluc_function(eta,cp2,cp3,R)

F=[];

cb = 0.5*erfc(eta);

a = (R-2-3*cb).*cp2 - (R-1)*cp3;

F = trapz(eta,a);

end

