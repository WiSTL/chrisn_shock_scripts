function [snr] = chris_snr_distribution(I,n)

[nx,ny] = size(I);

for i=1:nx-n
    for j=1:ny-n
        snr(i,j) = chris_snr(I(i:i+n,j:j+n),I(i:i+n,j:j+n));
    end
end

end