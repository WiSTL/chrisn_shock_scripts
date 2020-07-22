function [u,w] = chris_FFT_Poisson_solver(A,z,x)

w = A;

[dwdz,~] = chris_gradient(imgaussfilt(A,2),abs(z(2)-z(1)),abs(z(2)-z(1)));

u = abs(x(2)-x(1))*cumtrapz(-1*dwdz);

end
