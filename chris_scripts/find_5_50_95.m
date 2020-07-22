function [i5,i50,i95] = find_5_50_95(C)

[ny,~] = size(C);

C = medfilt2(C,[8 8]);

Cav = mean(C');

i5 = find(Cav<0.05);
if numel(i5)==0
i5 = 1;
else
i5 = i5(1);
end

i50 = find(Cav>0.49);
i50=i50(1);

i95 = find(Cav>0.95);
if numel(i95)==0
i95 = ny;
else
i95=i95(end);
end



end

