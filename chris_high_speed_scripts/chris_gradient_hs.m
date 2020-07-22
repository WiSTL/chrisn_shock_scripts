function [dudy,dudx,dudz] = chris_gradient_hs(u,dx,dy,dz)

[ny,nx,nz] = size(u);

for n = 1:nz
    [dudy(:,:,n),dudx(:,:,n)] = chris_gradient(squeeze(u(:,:,n)),dx,dy);
end

for n = 1:nx
    [~,dudz(:,n,:)] = chris_gradient(squeeze(u(:,n,:)),dz,dy);
end

% [nxx,nxy] = size(x);
% 
% if nxx~=nx
%    dx = x;
%    dy = y;
%    dz = z;
%    
%    xv = [1:nx]*dx;
%    yv = [1:ny]*dy;
%    zv = [1:nz]*dz;
%    
%    [y,x,z] = meshgrid(yv,xv,zv);
% end
% 
% %%% x derivative
% dudx=u;
% 
% for i = 3:nx-2
%     dudx(i,:,:) = (-u(i+2,:,:)+8*u(i+1,:,:)-8*u(i-1,:,:)+u(i-2,:,:))./(6*(y(i+1,:,:)-y(i-1,:,:)));
% end
% 
% for i = 1:2
%     dudx(i,:,:) = (u(i+1,:,:)-u(i,:,:))./(y(i+1,:,:)-y(i,:,:));
% end
% 
% for i = nx-1:nx
%     dudx(i,:,:) = (u(i-1,:,:)-u(i,:,:))./(y(i-1,:,:)-y(i,:,:));
% end
% 
% %%% y derivative
% dudy=u;
% 
% for i = 3:ny-2
%     dudy(:,i,:) = (-u(:,i+2,:)+8*u(:,i+1,:)-8*u(:,i-1,:)+u(:,i-2,:))./(6*(x(:,i+1,:)-x(:,i-1,:)));
% end
% 
% for i = 1:2
%     dudy(:,i,:) = (u(:,i+1,:)-u(:,i,:))./(x(:,i+1,:)-x(:,i,:));
% end
% 
% for i = ny-1:ny
%     dudy(:,i,:) = (u(:,i-1,:)-u(:,i,:))./(x(:,i-1,:)-x(:,i,:));
% end
% 
% %%% z derivative
% dudz=u;
% 
% for i = 3:nz-2
%     dudz(:,:,i) = (-u(:,:,i+2)+8*u(:,:,i+1)-8*u(:,:,i-1)+u(:,:,i-2))./(6*(z(:,:,i+1)-z(:,:,i-1)));
% end
% 
% for i = 1:2
%     dudz(:,:,i) = (u(:,:,i+1)-u(:,:,i))./(z(:,:,i+1)-z(:,:,i));
% end
% 
% for i = nz-1:nz
%     dudz(:,:,i) = (u(:,:,i-1)-u(:,:,i))./(z(:,:,i-1)-z(:,:,i));
% end
% 
% 
% % dudx(isnan(dudx))=0;
% % dudy(isnan(dudy))=0;
% % dudz(isnan(dudz))=0;


end

