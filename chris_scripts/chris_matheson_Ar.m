function [Q] = chris_matheson_Ar(x)

% x-rotameter reading
% Q-m3/s flow rate

Q = (2.475*10^(-9))*x^2 + (1.76*10^(-6))*x - 3.839*10^(-6);


end

