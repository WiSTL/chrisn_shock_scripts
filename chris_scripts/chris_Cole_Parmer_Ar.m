function [Q] = chris_Cole_Parmer_Ar(x)

% x-rotameter reading
% Q-m3/s flow rate

Q = (1.66667e-5) * (-0.002015*x^2 + (0.4208)*x - 1.895);


end

