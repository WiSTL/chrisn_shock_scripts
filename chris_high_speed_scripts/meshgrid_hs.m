function [E,Z,X,T] = meshgrid_hs(e,x,t,delta)

[ne,nt] = size(e);
nx = length(x);

for i = 1:ne
    for j = 1:nx
        for k = 1:nt
            E(i,j,k) = e(i,k);
            Z(i,j,k) = e(i,k)*delta(k);
            X(i,j,k) = x(j);
            T(i,j,k) = t(k);
        end
    end
end


end

