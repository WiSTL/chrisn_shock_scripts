function [u_new,v_new] = chris_PIV_remove_outliers(u,v,max,min)

[nx,ny,nimg] = size(u);

u(u>max) = NaN;
u(u<min) = NaN;

v(v>max) = NaN;
v(v<min) = NaN;

for n = 1:nimg
    for i = 2:nx-1
        for j = 2:ny-1
            if u(i,j,n)> 3*std([u(i,j,n),u(i+1,j,n),u(i-1,j,n),u(i,j+1,n),u(i,j-1,n),u(i+1,j+1,n),u(i-1,j+1,n),u(i+1,j-1,n),u(i-1,j-1,n)])
                u(i,j,n) = NaN;
            end
            if v(i,j,n)> 3*std([v(i,j,n),v(i+1,j,n),v(i-1,j,n),v(i,j+1,n),v(i,j-1,n),v(i+1,j+1,n),v(i-1,j+1,n),v(i+1,j-1,n),v(i-1,j-1,n)])
                v(i,j,n) = NaN;
            end
        end
    end 
    
    [u_new(:,:,n),~] = naninterp(u(:,:,n),v(:,:,n));
    [v_new(:,:,n),~] = naninterp(v(:,:,n),u(:,:,n));
    
end
end

