function [Sx,Sy,rx,ry] = chris_structure_fn(u,n)

    [nx, ny, nimg] = size(u);
    
    rx = [1:floor(nx/2)];
    
    Sx = rx*0;
    

    for i = 1:floor(nx/2)
        for j = rx
            Sx(j) = Sx(j) + mean(mean((abs(u(i+j,:,:) - u(i,:,:))).^n));
        end
    end
        
     Sx = Sx/floor(nx/2);
     
    ry = [1:floor(ny/2)];
    
    Sy = ry*0;
    

    for i = 1:floor(ny/2)
        for j = ry
            Sy(j) = Sy(j) + mean(mean((abs(u(:,i+j,:) - u(:,i,:))).^n));
        end
    end
        
     Sy = Sy/floor(ny/2);

end