function [F] = chris_zhou_mixed_mass_function(R)

eta = -50:0.05:50;

F=[];
[n] = length(R);

cb = 0.5*erfc(eta);

for i=1:n
    a = cb.*(1-cb).*(1+(R(i)-1)*cb);
    F(i) = 0.05*trapz(a);
end

end

