function [sig2,r] = chris_2D_delta_variance(u)

    [nx, ny] = size(u);
    um = nanmean(u(:));
    Sn = zeros(floor(sqrt((nx)^2+(ny)^2)),floor(sqrt((nx-1)^2+(ny-1)^2)));
    a = zeros(floor(sqrt((nx)^2+(ny)^2)),floor(sqrt((nx-1)^2+(ny-1)^2)));
    for i = 1:10:nx
        for j = 1:10:ny
            for k = 1:10:nx
               for l = 1:10:ny
                   if (i ~= k && j ~= l)
                  r = floor(sqrt((k-i)^2 + (l-j)^2));
                  x = floor(sqrt((i)^2 + (j)^2));
                  Sn(x,r) = Sn(x,r) +  u(k,l)-um;
                  a(x,r) = a(x,r) + 1;
                   end
               end
            end
        end
    end
    
    S = Sn./a;
    
    r = 1:sqrt((nx-1)^2+(ny-1)^2);
    dr = r(100)-r(99);
    
    sig = zeros(floor(sqrt((nx)^2+(ny)^2)),floor(sqrt((nx-1)^2+(ny-1)^2)));
    
    for L = 1:floor((2/3)*sqrt((nx-1)^2+(ny-1)^2))
        for i=1:floor(3*L/2)
            if i<L/2
                O = pi*((L/2)^-2)*1;
            else
                O = pi*((L/2)^-2)*(-0.125);
            end
            sig(:,L) = sig(:,L) + dr*(Sn(:,i)*O).^2;
        end
    end
    
    sig2 = nanmean(sig);
    
    r = r/1000;

end