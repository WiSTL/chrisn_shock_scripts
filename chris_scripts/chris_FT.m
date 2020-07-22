function [fh] = chris_FT(f,dx)

[nz,nx] = size(f);

k = [0:0.1:nx/2]*dx;
x = [1:nx]*dx;
z = [1:nz]*dx;

[Z,X] = meshgrid(z,x);

n = 0;

for i = k
    n=n+1;
    for m = 1:nz
        fh(m,n) = trapz(x,f(m,:).*exp(-2*pi*1j*i*x));
    end
end

end

