function [ia] = C_Contour(C,a)

[nx,ny] = size(C);

ia = [];

for i = 1:ny
    ias = find(C(:,i)<a);
    if isempty(ias)
        ias(i) = nx;
    else
        ia(i) = ias(1);
    end
end
end

