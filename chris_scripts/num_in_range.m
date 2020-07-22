function [n,u] = num_in_range(u,min,max)

u=u(:);
u = u(u>=min);
u = u(u<=max);

[n,~] = size(u);
end

