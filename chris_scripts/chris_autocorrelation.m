function [Rij,rx,Rji,ry] = chris_autocorrelation(u,v)


[nx, ny] = size(u);
    
rx = [1:floor(nx/2)];
    
Rij = rx*0;

 for i = 1:floor(nx/2)
        for j = rx
            Rij(j) = Rij(j) + mean(mean((u(i+j,:).*v(i,:))));
        end
 end
 
Rij = Rij/floor(nx/2);
Rij = Rij/Rij(1);

ry = [1:floor(ny/2)];
    
Rji = ry*0;

 for i = 1:floor(ny/2)
        for j = ry
            Rji(j) = Rji(j) + mean(mean((u(:,i+j).*v(:,i))));
        end
 end
 
Rji = Rji/floor(nx/2);
Rji = Rji/Rji(1);

end

