function [dudx,dudy] = chris_lele_gradient(u,dx,dy)

[nx,ny] = size(u);

for i = 1:nx
   dudy(i,:) = LeleD1_6(squeeze(u(i,:)),dy);
end

for i = 1:ny
   dudx(:,i) = LeleD1_6(squeeze(u(:,i)),dx); 
end

dudx(isnan(dudx))=0;
dudy(isnan(dudy))=0;

end

