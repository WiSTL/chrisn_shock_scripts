function [rx,Sp] = chris_structure_fn_p(u,p)

    [nx, ny] = size(u);
    
    if nx<ny
        rx = [1:floor(nx/2)];
    else
        rx = [1:floor(ny/2)];
    end
    
    Sp = rx*0;
    

    for i = rx
        for j = rx
            Sp(j) = Sp(j) + mean(mean((abs(u(i+j,:) - u(i,:))).^p));
        end
    end

    for i = rx
        for j = rx
            Sp(j) = Sp(j) + mean(mean((abs(u(:,i+j) - u(:,i))).^p));
        end
    end
        
     Sp = Sp/(floor(ny/2)+floor(ny/2));

end