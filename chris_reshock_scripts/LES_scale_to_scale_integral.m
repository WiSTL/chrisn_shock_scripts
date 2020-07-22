function [C,C2,r] = LES_scale_to_scale_integral(A,B,l,dl,dpxdx)

r = dl:dl:(l-dl);

r=r.^2;

n = 0;

C = [];

for i = r
    n = n+1;
    
    C(:,:,n) = chris_LES_filter_crop(chris_LES_filter_crop(A,sqrt(i)).*chris_LES_filter_crop(B,sqrt(i)),sqrt(l^2 - i)) - chris_LES_filter_crop(chris_LES_filter_crop(A,sqrt(i)),sqrt(l^2 - i)).*chris_LES_filter_crop(chris_LES_filter_crop(B,sqrt(i)),sqrt(l^2 - i));
    
end

C2 = squeeze(trapz(r/dpxdx^2,C,3));

end

