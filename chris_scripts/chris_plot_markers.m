function chris_plot_markers(x,y,z)

%%%%%%%%%%%%%%%%%%%%%
%
%   Creates marker plot of F=F(x,y) where the size and colour of 
%   the marker at the point (x,y) denotes the value of F
%
%%%%%%%%%%%%%%%%%%%%%

dx = nanmean(abs(diff(x)));
dy = nanmean(abs(diff(y)));

[ny nx] = size(z);

z(isnan(z))=0;
z(isinf(z))=0;

x(isnan(x))=0;
x(isinf(x))=0;
y(isnan(y))=0;
y(isinf(y))=0;

maxz = max(z(:));
minz = min(z(:));

hold on
grid on

for i = 1:nx
    for j = 1:ny
        if z(j,i)<=0
            d = 1 - (abs(minz) - abs(z(j,i)))/abs(minz);
            color = [0 0 1 0.9*d];
        else
             d = 1 - (abs(maxz) - abs(z(j,i)))/abs(maxz);
             color = [1 0 0 0.9*d];
        end
        xy = [x(i)-0.5*dx*d y(j)-0.5*dy*d d*dx d*dy];
        xy(isnan(xy))=0;
        xy(isinf(xy))=0;
        color(isnan(color))=0;
        rectangle('position',xy, 'FaceColor',color,'EdgeColor', 'none')
    end
end

xlim([min(x)-0.5 max(x)+0.5])
ylim([min(y)-0.5 max(y)+0.5])

hold off



end