function [C2] = shift_C(C)
[h,i5,i95] = concentration_centerline_5_95(C);
h = smooth(h);


[ny, nx] = size(C);

C2 = zeros(ny,nx);

for i = 1:nx    
   C2(1:ny-floor(h(i)-min(h)),i) = C(1+floor(h(i)-min(h)):end,i); 
end

end