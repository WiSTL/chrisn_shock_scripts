function [f,x,cp,a] = chris_hs_moment_func(C)

[nz,nx,nt] = size(C);

x = -10:0.05:10;

f = zeros(length(x),nz,nt);

for n = 1:nt
    
    Cm(:,n) = mean(squeeze(C(:,:,n)),2);
    
    [~,~,~,delta(n),~] = chris_C(squeeze(C(:,:,n)));
        
    for i = 1:nx
        Cp(:,i,n) = C(:,i,n) - Cm(:,n);
    end
    for j = 1:20
        cp(:,j,n) = mean(squeeze(Cp(:,:,n).^j),2);
        
        for k = 1:nz
            f(:,k,n) = f(:,k,n) + (1/factorial(j))*cp(k,j,n)*x.^(j-1)';
        end
    end  
    
    a(:,n) = trapz(cp(:,:,n))./delta(n);
end




end

