function [g,gamma,Rten,Ricci,Rs] = Surf_Find_Metric(A)

[~,Am,~,~] = C_fluctuation_595(A);
A = medfilt2(Am,[25 25]);

[nx,ny] = size(A);

[fx,fy] = gradient(A);

% fx = medfilt2(fx,[1 1]);
% fy = medfilt2(fy,[1 1]);

gxx = 1+(fx).^2;
gyy = 1+(fy).^2;
gzz = gxx./gxx;
gxy = fx.*fy;
gyx = gxy;
gxz = fx;
gzx = gxz;
gyz = fy;
gzy = gyz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%     Metric Tensor g ij
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = zeros(nx,ny,3,3);
gamma = zeros(nx-1,ny-1,3,3,3);
Rten = zeros(nx-2,ny-2,3,3,3,3);
Ricci = g;
Rs = 0*A(1:end-2,1:end-2);

g(:,:,1,1) = gxx;
g(:,:,1,2) = gxy;
g(:,:,2,1) = gxy;
g(:,:,2,2) = gyy;
g(:,:,2,3) = gyz;
g(:,:,3,2) = gyz;
g(:,:,3,3) = gzz;
g(:,:,1,3) = gxz;
g(:,:,3,1) = gzx;
g(:,:,3,2) = gzy;
g(:,:,2,3) = gzy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%     Christoffel Symbols gamma ijk
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:2
    for j=1:2
        for k=1:2
            a = diff(g(:,:,i,j),1,k);
            b = diff(g(:,:,i,k),1,j);
            c = diff(g(:,:,j,k),1,i);
            gamma(:,:,i,j,k) = a(1:nx-1,1:ny-1) + b(1:nx-1,1:ny-1) - c(1:nx-1,1:ny-1);
        end
    end
end

for i=1:2
    for j=1:2
        k=3;
        a = 0*gxx;
        b = diff(g(:,:,i,k),1,j);
        c = diff(g(:,:,j,k),1,i);
        gamma(:,:,i,j,k) = a(1:nx-1,1:ny-1) + b(1:nx-1,1:ny-1) - c(1:nx-1,1:ny-1);
    end
end

for i=1:2
    for k=1:2
        j=3;
        a = diff(g(:,:,i,j),1,k);
        b = 0*gxx;
        c = diff(g(:,:,j,k),1,i);
        gamma(:,:,i,j,k) = a(1:nx-1,1:ny-1) + b(1:nx-1,1:ny-1) - c(1:nx-1,1:ny-1);
    end
end

for j=1:2
    for k=1:2
        i=3;
        a = diff(g(:,:,i,j),1,k);
        b = diff(g(:,:,i,k),1,j);
        c = 0*gxx;
        gamma(:,:,i,j,k) = a(1:nx-1,1:ny-1) + b(1:nx-1,1:ny-1) - c(1:nx-1,1:ny-1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%     Riemann Curvature Tensor R ijkl
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:3
    for j=1:3
        for k=1:2
            for l=1:2
            a = diff(gamma(:,:,i,l,j),1,k);
            b = diff(gamma(:,:,i,k,j),1,l);
            c = 0*gamma(:,:,i,k,j);
            d=c;
            for nu = 1:3
                c = c+gamma(:,:,i,k,nu).*gamma(:,:,nu,l,j);
                d = d+gamma(:,:,i,l,nu).*gamma(:,:,nu,k,j);
            end  
            Rten(:,:,i,j,k,l) = a(1:nx-2,1:ny-2) - b(1:nx-2,1:ny-2) + c(1:nx-2,1:ny-2) - d(1:nx-2,1:ny-2);
            end
        end
    end
end

for i=1:2
    for j=1:2
        for k=1:2
            l=3;
            a = diff(gamma(:,:,i,l,j),1,k);
            b = 0*a;
            c = 0*gamma(:,:,i,k,j);
            d=c;
            for nu = 1:3
                c = c+gamma(:,:,i,k,nu).*gamma(:,:,nu,l,j);
                d = d+gamma(:,:,i,l,nu).*gamma(:,:,nu,k,j);
            end  
            Rten(:,:,i,j,k,l) = a(1:nx-2,1:ny-2) - b(1:nx-2,1:ny-2) + c(1:nx-2,1:ny-2) - d(1:nx-2,1:ny-2);
          
        end
    end
end

for i=1:2
    for j=1:2
        for l=1:2
            k=3;
            b = diff(gamma(:,:,i,k,j),1,l);
            a = 0*b;
            c = 0*gamma(:,:,i,k,j);
            d=c;
            for nu = 1:3
                c = c+gamma(:,:,i,k,nu).*gamma(:,:,nu,l,j);
                d = d+gamma(:,:,i,l,nu).*gamma(:,:,nu,k,j);
            end  
            Rten(:,:,i,j,k,l) = a(1:nx-2,1:ny-2) - b(1:nx-2,1:ny-2) + c(1:nx-2,1:ny-2) - d(1:nx-2,1:ny-2);
          
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%     Ricci Curvature Tensor R ij
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ricci = 0*g(1:nx-2,1:ny-2,:,:);

for i=1:3
    for j=1:3
        for p = 1:3
            Ricci(:,:,i,j) = Ricci(:,:,i,j)+Rten(:,:,p,i,p,j);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%     Ricci Scalar R
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rs = zeros(nx-2,ny-2);

for i = 1:3
    for j=1:3
        Rs(:,:) =Rs(:,:) + Ricci(:,:,i,j).*g(1:nx-2,1:ny-2,i,j);
    end
end

end
