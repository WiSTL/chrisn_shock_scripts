function [S, Smax, R, i5, i95] = C_shannon_entropy_595(C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cav = fliplr(mean(C'));

i5 = find(Cav<0.05);
if isempty(i5)
    i5 = 1;
else
    i5 = i5(end);
end
i5 = i5(end);

i95 = find(Cav>0.95);
if isempty(i95)
    i95 = length(Cav);
else
    i95 = i95(1);
end

C595 = C(i5:i95,1:end);

[nx,ny] = size(C595);

C5952 = 1-C595;

C595(C595<eps)=eps;
C5952(C5952<eps)=eps;

S=-sum(sum(C595.*log(C595) + C5952.*log(C5952) ))/(nx*ny)/log(2);

Im1 = sum(sum(C595))/(nx*ny);
Im2 = sum(sum(1-C595))/(nx*ny);

Smax = -(Im1*log(Im1)+Im2*log(Im2))/log(2);


R = sqrt(sum(sum(dot(gradient(C595.^0.5),gradient(C595.^0.5)))));


end

