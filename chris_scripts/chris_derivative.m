function [dudx] = chris_derivative(F,dx)

[dudx] = chris_lele_6(F,dx);

% [nx] = length(F);
% 
% if length(dx)==1
%    dx = ones(1,nx)*dx; 
% end
% 
% %%% x derivative
% dudx=F;
% 
% if nx>1
% 
% for i = 3:nx-2
%     dudx(i) = (-F(i+2)+8*F(i+1)-8*F(i-1)+F(i-2))/(12*dx(i));
% end
% 
% for i = 1:2
%     dudx(i) = (F(i+1)-F(i))/(dx(i));
% end
% 
% for i = nx-1:nx
%     dudx(i) = (F(i)-F(i-1))/(dx(i));
% end
% 
% else
%     dudx = 0;
% end

end

