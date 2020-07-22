function [h,F,f,k] = C_mixing_thickness_profile(C,K,pxcm)

h = 4*sum((C.*(1-C)));

[~,nx] = size(C);

k = (0:nx-1)/pxcm/nx;

f = (1/nx/pxcm)*smooth(fft(smooth(h)));

F = interp1(k,f,K);

end

