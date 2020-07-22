function [] = chris_scatter_plot(x,y,figno,num,holdon,xl,yl)

f = figure(figno);

set(gca, 'TickLabelInterpreter','latex'); 

if holdon
    hold on
else
    hold off
end

marker = ['+','s','x','v'];
colour = ['k','b','g','r'];

plot(x,y,marker(num),'MarkerSize',10,'MarkerEdgeColor',colour(num))

xlabel(xl);
ylabel(yl);

set(gca,'fontSize',18);

end

