function [h,F,f,k] = C_scalar_power_spectra(C,K,pxcm)

h = C;

[i5,~,i95] = find_5_50_95(C);

[~,nx] = size(C);

[f,k] = power_spectra_1D(h(i5:i95,:),1,0,1);

f = mean(f);
k = k/pxcm;

F = (interp1(k,f,K));

end

