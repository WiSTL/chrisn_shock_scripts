function [dudy,dudx,dudz] = chris_gradient_hs_uniform(u,dx,dy,dz)

[~,nx,nz] = size(u);

for k = 1:nz
    [dudy1(:,:,k),dudx(:,:,k)] = chris_gradient(squeeze(u(:,:,k)),dx,dy);
end

for i = 1:nx
    [dudy2(:,i,:),dudz(:,i,:)] = chris_gradient(squeeze(u(:,i,:)),dz,dy);
end

dudy = 0.5*(dudy1+dudy2);

end

