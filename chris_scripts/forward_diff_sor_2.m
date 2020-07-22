function [p,itercount,difference] = forward_diff_sor_2(a,b,c,p,r,nx,ny)
%-----------------------------------------------------------------------
%SOR solver
%-----------------------------------------------------------------------
%initialise
err = 1.0e-4; 
difference = 1.0; 
itermax = 10000; 
itercount = 0;

b = medfilt2(b);
bi = 1./b;

bi(isnan(bi)) = 0;
bi(isinf(bi)) = 0;

a = medfilt2(a);
c = medfilt2(c);
r = medfilt2(r);
bi = medfilt2(bi);

figure(1)
imagesc(a)
colorbar

figure(2)
imagesc(c)
colorbar

figure(3)
imagesc(r)
colorbar

figure(4)
imagesc(b)
colorbar

figure(5)
imagesc(bi)
colorbar

while ( difference > err && itercount < itermax )

itercount = itercount + 1;

if( itercount > itermax ) 
    disp( 'method did not converge' ); 
end

pold = p;

%sor method 
difference = 0; 
for j = 2:ny-1 
    for i = 2:nx-1
        
        if abs(b(i-1,j))<10 || abs(b(i,j))<10 || abs(b(i,j-1))<10
            p(i,j) = p(i,j);
        else
            p(i,j) = r(i,j).*bi(i,j) + a(i,j)*p(i,j-1).*bi(i,j) + c(i,j)*p(i-1,j).*bi(i,j); 
            p(isnan(p))=0;
            p(isinf(p))=0;
        end
         difference = difference + abs(p(i,j)-pold(i,j));
    end
end
%normalise
clc
difference = difference/((nx-1)*(ny-1))

figure(6)
imagesc(p)
colorbar
pause(0.00001)
end
end