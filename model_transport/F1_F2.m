function [F1,F2] = F1_F2(n,eta_L)

F1 = 1.74*(n.^-6).*(eta_L.^-1.5) +1.178;

F2 = 131.1 - 12.01*n - 871.5*eta_L;

end

