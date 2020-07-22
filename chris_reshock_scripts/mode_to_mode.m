[uh,~,~] = fft2d(u);
[vh,kx,ky] = fft2d(v);

[Kx,Ky] = meshgrid(kx,ky);

Suu111 = imag((Kx(1,1)*uh(1,2)+Ky(1,1)*vh(1,2))*(uh(1,1)^2+ vh(1,1)^2));

[ny,nx] = size(uh);

T = zeros(nx,nx);

no = 0;

for j = 1:10:ny
    no = no+1;
    for n = 1:nx
        for m = 1:nx
            ki = n;
            pi = m;
            qi = 1+mod(n+m,nx);
        
                T(m,n,no) = imag((Kx(j,ki)*uh(j,qi)+Ky(j,ki)*vh(j,qi))*(uh(j,ki)*uh(j,pi)+ vh(j,ki)*vh(j,pi)));
        end
    end
end