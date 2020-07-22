function [P,x,y,cov,r,P_xy,P_yx] = chris_jpdf_moments(U,V,n,MIN_XY,MAX_XY)



[P,x,y] = chris_jpdf(U,V,n,MIN_XY,MAX_XY);

[ny,nx] = size(P);

x1 = x(1,:);
y1 = y(:,1);

%%%
% Independent Probability Density Functions
%%%
Pu = trapz(y1,P);
Pv = trapz(x1,P');

%%%
% Means
%%%
meanu = trapz(x1,x1.*Pu);
meanv = trapz(y1,y1.*Pv');

%%%
% Variances
%%%
varu = trapz(x1,((x1-meanu).^2).*Pu);
varv = trapz(y1,((y1-meanv).^2).*Pv');

%%%
% Covariance
%%%
cov = trapz(x1,trapz(y1,x.*y.*P)) - meanu*meanv;

%%%
% correlation
%%%
r = cov/(sqrt(varu)*sqrt(varv));

%%%
% Conditional PDFs
%%%
for i = 1:nx
   P_xy(:,i) = P(:,i)./Pu(i); 
end

for i = 1:ny
   P_yx(i,:) = P(i,:)./Pv(i); 
end

end

