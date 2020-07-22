function [Sx,rx,dA,Am] = chris_no_structure_fn(A,n,dpxdx)

    [nz,nx,~] = size(A);
    
    if mod(nx,2)==1
        nx=nx-1;
    end
    
    rx = [1:floor(nx/2)];
    
    Sx = 0*A(:,1:ceil(nx/2));
    
    for i = 1:ceil(nx/2)
        for j = rx
            if (i+j)<=nx
                Sx(:,j) = Sx(:,j) + (abs(A(:,i+j) - A(:,i))).^n;
                dA(:,i,j) = A(:,i+j) - A(:,i);
                Am(:,i,j) = A(:,i+j) + A(:,i);
            else
                Sx(:,j) = Sx(:,j) + (abs(A(:,2*nx-(i+j)+1) - A(:,i))).^n; 
                dA(:,i,j) = A(:,i+j) - A(:,i);
                Am(:,i,j) = A(:,i+j) + A(:,i);
            end            
        end
    end
    
     Sx = Sx/floor(nx/2);
     rx = (rx-1)/dpxdx;
     
%      for i = 1:ceil(nx/2)
%         [dAdz(:,i,:),dAdr(:,i,:)] = chris_gradient(squeeze(dA(:,i,:)),1/dpxdx,0.01);
%      end
%      
%      for i = 1:nz
%         [dAdx(i,:,:),~] = chris_gradient(squeeze(dA(i,:,:)),1/dpxdx,1/dpxdx);
%      end
%     
%      xzn = n*(n-1)*squeeze(nanmean(((dA).^(n-2)).*(dAdz.^2),2));
%      xrn = n*(n-1)*squeeze(nanmean(((dA).^(n-2)).*((dAdr.^2)+(dAdx.^2)),2));
%      
end