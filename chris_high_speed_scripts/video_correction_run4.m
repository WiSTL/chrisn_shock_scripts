
C_corrected = [];
c1=[];
c2=[];

ak=[];
al=[];

% for j=1:61
%    top(j) = 50; 
%    bottom(j) = floor(700 - (j-1)*(700-200)/(61-1)); 
% end
% 
% top(62:90)=50;
% bottom(62:90) = 200;
% 
% for j=91:147
%    top(j) = 50; 
%    bottom(j) = floor(200 - (j-91)*(200-300)/(147-91)); 
% end

bottom = 500;

for j = 1:147
    c1(:,:,j) = lasertransform(squeeze(c(:,:,j)), origin, 'forward');
end


for j = 1:147
    for i = 1:278
        af = squeeze(c1(:,i,j));
        faf(:,i,j) = filter(w/sum(w),1,af);
%         [~,bottom] = max(smooth(squeeze(faf(:,i,j)),300));
        S0 = squeeze(faf(bottom,i,j));
        x=[top(j):bottom]'; %rowlow]';%
        y=smooth(squeeze(faf(top(j):bottom,i,j)),10);%rowlow);%
        logy=real(log(y));
    
        infIndexes = isinf(logy);
        logy(infIndexes) = [];
        x = x(~infIndexes);
    
        linefit=polyfit(x(:),logy(:),1);
        intercept=real(exp(linefit(2)));
        abs_col=abs(real(linefit(1)));
        
        C_corrected(:,i,j) = squeeze(faf(:,i,j))./(S0 + abs_col*cumtrapz(squeeze(faf(:,i,j))));
        
    end
end

for j = 1:147
    c2(:,:,j) = lasertransform(squeeze(C_corrected(:,:,j)), origin, 'inverse');
end


% for i = 1:97
%     i
%     C_corrected(:,:,i) = chris_c_correction(squeeze(c(:,:,i)),[top(i) bottom(i) 1 278], origin);
%     
% end