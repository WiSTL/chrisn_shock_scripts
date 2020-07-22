function [x,y1,y2,y3] = chris_plot_mean_std(x,ym,c,n)

[ny,nx] = size(ym);

y1 = smooth(nanmean(ym'),n)';
y2 = smooth(nanmean(ym')+2*nanstd(ym')/sqrt(nx),n)';
y3 = smooth(nanmean(ym')-2*nanstd(ym')/sqrt(nx),n)';

X=[x,fliplr(x)];                %#create continuous x value array for plotting
Y=[y2,fliplr(y3)];              %#create y values for out and then back

h = fill(X,Y,c); 

set(h,'FaceAlpha',0.3);
set(h,'edgealpha',0)

hold on;
plot(x,y1,c,'linewidth',2);
hold off


end

