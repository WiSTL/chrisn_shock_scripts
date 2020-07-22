function [p,itercount,difference] = forward_diff_sor_4(a,b,p,r,nx,ny,dx,dz)
%-----------------------------------------------------------------------
%SOR solver
%-----------------------------------------------------------------------
%initialise
err = 1.0e-4; 
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

a(isnan(a))=0;
b(isnan(b))=0;
r(isnan(r))=0;

while ( difference > err && itercount < itermax )

itercount = itercount + 1;

if( itercount > itermax ) 
    disp( 'method did not converge' ); 
end

p(isnan(p))=0;

pold = p;

ri = r./a;
bi = b./a;

ri(isnan(ri))=0;
bi(isnan(bi))=0;

ri(isinf(ri))=0;
bi(isinf(bi))=0;

%sor method 
difference = 0.0; 
for n = 4:nx-4
    for m = 4:ny-4 
        if abs(a(n,m))< 10
            p(n+1,m) = (-1/39)*(-39*pold(n-1,m) + 12*(pold(n+2,m) - pold(n-2,m)) - 5*(pold(n+3,m) - pold(n-3,m)));
        else
            p(n+1,m) = (-1/39)*(-39*pold(n-1,m) + 12*(pold(n+2,m) - pold(n-2,m)) - 5*(pold(n+3,m) - pold(n-3,m))) - (bi(n,m)/39)*(dx/dz)*(39*(pold(n,m+1) - pold(n,m-1)) + 12*(pold(n,m+2) - pold(n,m-2)) - 5*(pold(n,m+3) - pold(n,m-3))) + ri(n,m)*(96*dx/39);
        end
        
        difference = difference + abs(p(n,m)-pold(n,m));
        
    end
end
%normalise
clc
difference = difference

figure(4)
imagesc(p)
colorbar
pause(0.1)
end
end