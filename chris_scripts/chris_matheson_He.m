function [Q] = chris_matheson_He(x)

% x-rotameter reading
% Q-m3/s flow rate

Q = (2.035*10^(-8))*x^2 + (3.28*10^(-6))*x - 1.518*10^(-5);


end

