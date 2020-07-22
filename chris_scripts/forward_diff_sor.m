function [p,itercount,difference] = forward_diff_sor(a,b,c,p,r,nx,ny)
%-----------------------------------------------------------------------
%SOR solver
%-----------------------------------------------------------------------
%initialise
err = 1.0e-4; 
difference = 1.0; 
itermax = 10000; 
itercount = 0;

b(abs(b)<100)=0;
bi = 1./b;

bi(isnan(bi)) = 0;
bi(isinf(bi)) = 0;

bi = medfilt2(bi,[20 20]);

a = a.*bi;
c = c.*bi;

r = r.*bi;

a(isnan(a)) = 0;
a(isinf(a)) = 0;

c(isnan(c)) = 0;
c(isinf(c)) = 0;

r(isnan(r)) = 0;
r(isinf(r)) = 0;

a = medfilt2(a,[10 10]);
c = medfilt2(c,[20 20]);
r = medfilt2(r);

figure(1)
imagesc(a)
colorbar

figure(2)
imagesc(c)
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
p2 = p;

%sor method 
difference = 0.0; 
for j = 2:ny-1 
    for i = 2:nx-1  
        p(i,j) = r(i,j) + (a(i,j)+a(i-1,j))*p(i-1,j) + (c(i,j)+c(i,j-1))*p(i,j-1); 
        p(isnan(p))=0;
        p(isinf(p))=0;
        p2(i,j) = r(i,j) + (a(i,j)+a(i+1,j))*p(i+1,j) + (c(i,j)+c(i,j+1))*p(i,j+1);
        difference = difference + abs(p2(i,j)-pold(i,j)); 
    end
end
%normalise
clc
difference = difference/((nx-1)*(ny-1))

p = p2;

figure(4)
imagesc(p)
colorbar
pause(0.00001)
end
end