function [dudy,dudx] = chris_gradient(u,x,y)

[ny,nx] = size(u);

for i = 1:ny
    dudx(i,:) = chris_derivative(squeeze(u(i,:)),x);
end

for i = 1:nx
    dudy(:,i) = chris_derivative(squeeze(u(:,i)),y);
end



% [nx,ny] = size(u);
% 
% [nxx,nxy] = size(x);
% 
% if nxx~=nx || nxy~=ny
%    dx = x;
%    dy = y;
%    
%    xv = [1:nx]*dx;
%    yv = [1:ny]*dy;
%    
%    [y,x] = meshgrid(yv,xv);
% end    
% 
% %%% x derivative
% dudy=u;
% if nx>4
% for i = 3:nx-2
%     dudy(i,:) = (-u(i+2,:)+8*u(i+1,:)-8*u(i-1,:)+u(i-2,:))./(6*(x(i+1,:)-x(i-1,:)));
% end
% end
% 
% for i = 1:2
%     dudy(i,:) = (u(i+1,:)-u(i,:))./(x(i+1,:)-x(i,:));
% end
% 
% for i = nx-1:nx
%     dudy(i,:) = (u(i-1,:)-u(i,:))./(x(i-1,:)-x(i,:));
% end
% 
% %%% y derivative
% dudx=u;
% if ny>4
% for i = 3:ny-2
%     dudx(:,i) = (-u(:,i+2)+8*u(:,i+1)-8*u(:,i-1)+u(:,i-2))./(6*(y(:,i+1)-y(:,i-1)));
% end
% end
% 
% for i = 1:2
%     dudx(:,i) = (u(:,i+1)-u(:,i))./(y(:,i+1)-y(:,i));
% end
% 
% for i = ny-1:ny
%     dudx(:,i) = (u(:,i-1)-u(:,i))./(y(:,i-1)-y(:,i));
% end
% 
% 
% dudx(isnan(dudx))=0;
% dudy(isnan(dudy))=0;
% 
% dudx = dudx;
% dudy = dudy;

end

