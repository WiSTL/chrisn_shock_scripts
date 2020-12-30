
nx = 20;
ny = 32;
nz = 32;

dx = 2*pi/(nx-1);
dy = 2*pi/(ny-1);
dz = 2*pi/(nz-1);

x = [0:dx:2*pi];
y = 0:dy:2*pi;
z = 0:dz:2*pi;

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