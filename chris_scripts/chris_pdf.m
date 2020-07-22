function [A,x] = chris_pdf(U,interval,minU,maxU)

U = U(:);

n = size(U);

if isnan(minU)
minU = min(U);
end

if isnan(maxU)
maxU = max(U);
end

x = minU:interval:maxU-interval;
A = zeros(1,length(x));


 U = round(U,3);
 x = round(x,3);

n=0;
for i=x
    n=n+1;
    [A(n),~] = num_in_range(U,i,i+interval);    
end

A = medfilt1(A,10);

dx = mean(diff(x));

normA = sum(A*dx);

A = A/normA;

end

