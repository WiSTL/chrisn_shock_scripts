function [T,numT] = chris_fourier_mode_energy_transfer(u,v)

[uh,~,~] = fft2d(u);
[vh,kx,ky] = fft2d(v);

[Kx,Ky] = meshgrid(kx,ky);

[ny,nx] = size(uh);

T = zeros(floor(sqrt(nx^2+ny^2)),floor(sqrt(nx^2+ny^2)));
numT = zeros(floor(sqrt(nx^2+ny^2)),floor(sqrt(nx^2+ny^2)));

for in = 1:5:nx
    for jn = 1:5:ny
       for im = 1:20:nx
           for jm = 1:20:ny
               n = floor(sqrt(in^2+jn^2));
               m = floor(sqrt(im^2+jm^2));
               
               q =1+mod(n+m,min(nx,ny));
               for jq = 1:q-1
                   iq = sqrt(q^2-jq^2);
                   if (iq == floor(iq))
                       numT(m,n) = numT(m,n)+1;
                       T(m,n) = T(m,n) + imag((Kx(jn,in)*uh(jq,iq)+Ky(jn,in)*vh(jq,iq))*(uh(jn,in)*uh(jm,im)+ vh(jn,in)*vh(jm,im)));
                   end
               end
           end
       end
    end
end
end