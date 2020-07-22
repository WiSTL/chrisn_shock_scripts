function [Cfrms,Am,Cm,Cym,h595] = C_fluctuation_595(A)

Aav = fliplr(mean(A'));

i5 = find(Aav<0.05);
if isempty(i5)
    i5 = 1;
else
    i5 = i5(end);
end
i5 = i5(end);

i95 = find(Aav>0.95);
if isempty(i95)
    i95 = length(Aav);
else
    i95 = i95(1);
end

h595 = abs(i5-i95);

A595 = A(i5:i95,1:end);

Cm = mean(mean(A595));

[nx ny] = size(A595);

meanA = mean(A595');

Cym = ones(nx,ny);

for i=1:ny
    Cym(:,i) = meanA';
end


Am = A595;

for i = 1:ny
   Am(:,i) = A595(:,i)- meanA';
end

Cfrms = sqrt(sum(sum(Am.^2))/(nx*ny));

end