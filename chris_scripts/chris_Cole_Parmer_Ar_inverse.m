function [x] = chris_Cole_Parmer_Ar_inverse(Q)

% x-rotameter reading
% Q-m3/s flow rate

Q = Q/(1.66667e-5);

x = 0.0999*Q^2 + 1.077*Q - 0.4097;


end

