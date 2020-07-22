function [x] = chris_matheson_He_inverse(Q)

% x-rotameter reading
% Q-m3/s flow rate

x = (-6.46*10^(7))*Q^2 + (2.122*10^(5))*Q + 7.569;


end