function [p,itercount,difference] = forward_diff_sor_3(a,b,p,r,nx,ny)
%-----------------------------------------------------------------------
%SOR solver
%-----------------------------------------------------------------------
%initialise
err = 1.0e-3; 
difference = 1.0; 
itermax = 10000; 
itercount = 0;

figure(1)
imagesc(a)
colorbar

figure(2)
imagesc(b)
colorbar

figure(3)
imagesc(r)
colorbar

while ( difference > err && itercount < itermax )

itercount = itercount + 1;

if( itercount > itermax ) 
    disp( 'method did not converge' ); 
end

pold = p;

%sor method 
difference = 0.0; 
for i = 2:nx-1
    for j = 2:ny-1 
        if abs(a(i,j))< 100
            p(i,j+1) = pold(i,j+1);
        else
            p(i,j+1) = pold(i,j-1) + r(i,j)./a(i,j) - (b(i,j)/a(i,j))*(pold(i+1,j)-pold(i-1,j));
        end
    end
end
for j = 2:ny-1
    for i = 2:nx-1  
        if abs(b(i,j))< 100
            p(i+1,j) = pold(i+1,j);
        else
            p(i+1,j) = pold(i-1,j) - r(i,j)./b(i,j) + (a(i,j)/b(i,j))*(pold(i,j+1)-pold(i,j-1));
        end
        difference = difference + abs(p(i,j)-pold(i,j));
    end
end
%normalise
clc
difference = difference/(nx*ny)

figure(4)
imagesc(p)
colorbar
pause(0.1)
end
end