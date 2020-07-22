function [h] = jet_velocity_centerline_05(V)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[nx,ny] = size(V);

h = zeros(1,ny);

% for i=1:ny
%     
%     i5 = find(smooth(C(:,i),200)<0.05);
%     if isempty(i5)
%     	i5 = 1;
%     else
%         i5 = i5(1);
%     end
%     i5 = i5(end);
% 
%     i95 = find(smooth(C(:,i),200)>0.95);
%     if isempty(i95)
%         i95 = length(C(:,i));
%     else
%         i95 = i95(end);
%     end
% 
% h(i) = abs(i95-i5);
% 
% end

for i=1:ny
    
    vhalf = max(V(:,i))/2;
    
    for j=1:nx
        if V(j,i)>vhalf
            i95 = j;
            break
        end
    end
    
    for j=nx:-1:1
        if V(j,i)>vhalf
            i5 = j;
            break
        end
    end
    
    h(i) = (i95+i5)/2;
    
end


end

