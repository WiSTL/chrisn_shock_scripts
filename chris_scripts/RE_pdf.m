function [pdf,l] = RE_pdf(C)

[nx,ny] = size(C);

l = 1:10:ny;

n=0;

for i = 1:10:ny
   n=n+1;
   pdf(n) = i*abs(max(max(imresize(C,[i i]))));
    
end
end

