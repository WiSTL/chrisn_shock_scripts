
nx = 1;
ny = 128;
nz = 256;

dx = 2*pi/(1024-1);
dy = 2*pi/(1024-1);
dz = 2*pi/(1024-1);

x = [0:dx:dx*(nx-1)];
y = 0:dy:dy*(ny-1);
z = 0:dz:dz*(nz-1);

rhom = [];
uxm = [];
uym = [];
uzm = [];

n = 0;

for i = 1:nx
    for j = 1:ny
        for k = 1:nz
           n=n+1;
            
           rhom(i,j,k) = rho(n);
           uxm(i,j,k) = ux(n); 
           uym(i,j,k) = uy(n); 
           uzm(i,j,k) = uz(n); 
        end
    end
end