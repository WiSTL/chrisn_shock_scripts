function [Sc1D,kx] = Chris_Ensemble_1D_Spectra(A,dx)


[~,~,nimg] = size(A);

[Sc1D,kx] = power_spectra_1D(A(:,:,1),dx);

for i = 2:nimg
    [Sc1Da,kx] = power_spectra_1D(A(:,:,i),dx);
    Sc1D = Sc1D+Sc1Da;
end

Sc1D = Sc1D/nimg;


end

