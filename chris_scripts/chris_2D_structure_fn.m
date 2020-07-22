function [S,r] = chris_2D_structure_fn(u,n)

    [nx, ny] = size(u);
    Sn = zeros(1,floor(1000*sqrt((nx-1)^2+(ny-1)^2)));
    a = zeros(1,floor(1000*sqrt((nx-1)^2+(ny-1)^2)));
    for i = 1:10:nx
        for j = 1:10:ny
            for k = 1:5:nx
               for l = 1:5:ny
                   if (i ~= k && j ~= l)
                  r = floor(1000*sqrt((k-i)^2 + (l-j)^2));
                  
                  Sn(r) = Sn(r) +  (abs(u(k,l) - u(i,j))).^n;
                  a(r) = a(r) + 1;
                   end
               end
            end
        end
    end
    S = Sn./a;
    r = 1:1000*sqrt((nx-1)^2+(ny-1)^2);
    r = r/1000;

    n = 0;
    for i = 1:length(S)
       if isnan(S(i))
       else
           n=n+1;
           S(n) = S(i);
           r(n) = r(i);
       end
    end
    
    S = S(1:n);
    r = r(1:n);
    
end