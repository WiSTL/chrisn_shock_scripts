function [Sux,Suy,Svx,Svy,Sxy,Syx,rux,ruy,rvx,rvy,rxy,ryx] = chris_2o_structure_fn(u,v)

    [nx, ny, nimg] = size(u);
    
    rux = [1:floor(nx/2)];
    
    Sux = rux*0;
    

    for i = 1:floor(nx/2)
        for j = rux
            Sux(j) = Sux(j) + mean(mean(((u(i+j,:,:) - u(i,:,:))).^2));
        end
    end
        
     Sux = Sux/floor(nx/2);
     
    ruy = [1:floor(ny/2)];
    
    Suy = ruy*0;
    

    for i = 1:floor(ny/2)
        for j = ruy
            Suy(j) = Suy(j) + mean(mean(((u(:,i+j,:) - u(:,i,:))).^2));
        end
    end
        
     Suy = Suy/floor(ny/2);
     
     rvx = [1:floor(nx/2)];
    
    Svx = rvx*0;
    

    for i = 1:floor(nx/2)
        for j = rvx
            Svx(j) = Svx(j) + mean(mean(((v(i+j,:,:) - v(i,:,:))).^2));
        end
    end
        
     Svx = Svx/floor(nx/2);
     
    rvy = [1:floor(ny/2)];
    
    Svy = rvy*0;
    

    for i = 1:floor(ny/2)
        for j = rvy
            Svy(j) = Svy(j) + mean(mean(((v(:,i+j,:) - v(:,i,:))).^2));
        end
    end
        
     Svy = Svy/floor(ny/2);
     
    rxy = [1:floor(ny/2)];
    
    Sxy = rxy*0;
    

    for i = 1:floor(ny/2)
        for j = rxy
            Sxy(j) = Sxy(j) + mean(mean(((u(:,i+j,:) - u(:,i,:))).*(v(:,i+j,:) - v(:,i,:))));
        end
    end
        
     Sxy = Sxy/floor(ny/2);

     ryx = [1:floor(nx/2)];
    
    Syx = ryx*0;
    

    for i = 1:floor(nx/2)
        for j = ryx
            Syx(j) = Syx(j) + mean(mean(((u(i+j,:,:) - u(i,:,:))).*(v(i+j,:,:) - v(i,:,:))));
        end
    end
        
     Syx = Syx/floor(nx/2);
     
end