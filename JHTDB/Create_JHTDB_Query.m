

nx = 20;
ny = 32;
nz = 32;

dx = 2*pi/(nx-1);
dy = 2*pi/(ny-1);
dz = 2*pi/(nz-1);

x = [0:dx:2*pi];
y = 0:dy:2*pi;
z = 0:dz:2*pi;

A = [];
B = [];
C = [];

n=0;

for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            n=n+1;
            
            A(n) = x(i);
            B(n) = y(j);
            C(n) = z(k);
            
        end
    end
end

fileID = fopen('query.txt','w');
formatSpec = '%1.10f %1.10f %1.10f \n';
fprintf(fileID,formatSpec,A,B,C);
fclose(fileID);

