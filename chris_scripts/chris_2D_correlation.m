function [C,r] = chris_2D_correlation(u,w)

    [nx, ny] = size(u);
    Sn = zeros(1,floor(1000*sqrt((nx-1)^2+(ny-1)^2)));
    a = zeros(1,floor(1000*sqrt((nx-1)^2+(ny-1)^2)));
    for i = 1:10:nx
        for j = 1:10:ny
            for k = 1:5:nx
               for l = 1:5:ny
                   if (i ~= k && j ~= l)
                  r = floor(1000*sqrt((k-i)^2 + (l-j)^2));
                  Sn(r) = Sn(r) +  ((u(k,l)*w(i,j)));
                  a(r) = a(r) + 1;
                   end
               end
            end
        end
    end
    C = Sn./a;
    r = 1:1000*sqrt((nx-1)^2+(ny-1)^2);
    r = r/1000;
    
    n = 0;
    for i = 1:length(C)
       if isnan(C(i))
       else
           n=n+1;
           C(n) = C(i);
           r(n) = r(i);
       end
    end
    
    C = C(1:n);
    r = r(1:n);

end