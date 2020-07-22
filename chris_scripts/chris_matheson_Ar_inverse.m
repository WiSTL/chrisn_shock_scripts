function [x] = chris_matheson_Ar_inverse(Q)

% x-rotameter reading
% Q-m3/s flow rate

x = (-2.39*10^(8))*Q^2 + (5.418*10^(5))*Q + 2.638;


end