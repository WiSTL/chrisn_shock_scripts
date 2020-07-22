function [Cnew] = chris_diffusion_correction(C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Corrects for the differential diffusion
%%% that gives values of C above 1 after
%%% correction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx,ny,nimg] = size(C);

for n = 1:nimg
    Ci = C(:,:,n);
    
   f = mean(Ci');
   h = mean(f(end-100:end));
   for i = 1:nx
       for j = 1:ny
          Ci(i,j) = Ci(i,j) - h;
       end
   end 
   f = mean(Ci(1:500,:));
   
   for i = 1:ny
       Ci(:,i) = Ci(:,i)/f(i);
   end
   
%    g = mean(f(1:50));
%    
%    Ci = Ci/g;
   Ci(Ci>1) = 1;
   
   Cnew(:,:,n) = Ci;
end

Cnew(Cnew>1)=1;
Cnew(Cnew<0)=0;

end
