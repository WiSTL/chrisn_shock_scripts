function [I1] = chris_binning(I,xbins,ybins)

[ny,nx] = size(I);

dx = floor(nx/xbins);
dy = floor(ny/ybins);


for i = 1:xbins
    for j=1:ybins
        I1(i,j) = nansum(nansum(I((j-1)*dy+1:j*dy,(i-1)*dx+1:i*dx)));
    end
end

end

