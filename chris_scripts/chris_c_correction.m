function [imgFinal,X,imgT,abs_col,intercept,logy,x] = chris_c_correction(img,params,origin)

rowtop=params(1);rowlow=params(2); % rows of uniform seeding, laser attenuation correction
colleft=params(3); colright=params(4); % columns of uniform seeding


Xlimits = [0.05, 0.95]; % Relative acetone concentration lim for img crop
pxcm = 63;        % Pixels per cm in imaging plane of input images

pad_px = 1024;

imgRaw = img;
imgRaw = imgRaw/max(imgRaw(:));

imgRaw = double(imgRaw);
img=imgRaw;
[nrow ncol]=size(imgRaw);
img(img<0) = 0;

%% locate origin of laser striations
% [origin, lines] = laserorigin_JGO(img(rowtop:rowlow,colleft:colright));
% hold on;
% title(['Verify that green lines are collinear with laser striations. '...
%        'Modify laserorigin.m if necessary.']);
% for k = 1:length(lines)
%    xy = [lines(k).point1; lines(k).point2];
%    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% end
% hold off
% pause(0.5);

% close
%% Transform image so laser striations are vertical
imgT = lasertransform(img, origin, 'forward');

[nrow, ncol]=size(imgT);

%% Obtain uniform region for given image
% figure(1)
% imagesc(imgT)
% caxis([0 1])
% [x, y] = ginput;
% 
% x_m = colleft:colright;
% 
% y2 = floor(interp1(x,y,x_m,'pchip'));
% 
% hold on
% plot(x_m,y2,'Linewidth',3);
% hold off
%% Correct for laser attenuation
imgCorr=imgT*0;
imgTf=imgaussfilt(imgT,3);
imgTf=imgT;
abs_col=zeros(ncol,1);
intercept=zeros(ncol,1);
for k=1:ncol
    intcolumn=imgTf(:,k);
    
    x=[rowtop:rowlow]';
    y=intcolumn(rowtop:rowlow);
    logy=real(log(y));
    
    infIndexes = isinf(logy);
    logy(infIndexes) = [];
    x = x(~infIndexes);
    
    linefit=polyfit(x,logy,1);
    intercept(k)=real(exp(linefit(2)));
    abs_col(k)=real(linefit(1));
%    imgCorr(:,k)=imgT(:,k)./intabs'; % divides by curvefit for correction
% this makes mixed regions lighter (in error) because it assumes entire
% region is uniformly seeded
end

% obtain correction by dividing by integral of signal (Weber thesis eq 4.8)
correct=zeros(nrow,1);

abscol = nanmean(abs_col);

for k=1:ncol
    correct=intercept(k)+abscol*cumtrapz(imgTf(:,k));
    imgCorr(:,k)=imgT(:,k)./correct; % corrects for attenuation
end
imgCorr(isnan(imgCorr)) = 0; 

%% Correct for non-uniform beam profile and artifacts from notch
imgTemp=imgCorr;
% imgTf=imgaussfilt(imgTemp,3);
imgTf=imgTemp;
% meanuniform=mean2(imgTemp(rowtop:rowlow,colleft:colright));
for k=1:ncol
    imgTemp(:,k)=imgTemp(:,k)/mean(imgTf(rowtop:rowlow,k));
end
imgCorr=imgTemp;
imgCorr(isnan(imgCorr)) = 0; 
%% Transform back to Cartesian frame
imgCorr(isnan(imgCorr))=0;
imgCorr(isinf(imgCorr))=0;
imgCorr = imgCorr/nanmean(nanmean(imgCorr(rowtop:rowlow, colleft:colright)));
imagesc(imgCorr);
imgFinal = lasertransform(imgCorr, origin, 'inverse');
% imgFinal = imgFinal/nanmean(nanmean(imgFinal(rowtop:rowlow, colleft:colright)));

X = acetoneX(medfilt2(imgCorr), abs_col, pxcm, Xlimits, pad_px);

end

