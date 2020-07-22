function  [p] = chris_plot_anim(x,y,range,dt,xl,yl,fn)

[ni,nt] = size(y);

figure
hold on

for i = 1:nt
    if i>1
       p(i-1). Color(4) = 0.75;
       p(i-1).LineWidth = 1.75;
       
%        q(i-1). Color(4) = 0.75;
%        q(i-1).LineWidth = 1.75;
    end
    if i>2
       p(i-2). Color(4) = 0.5;
       p(i-2).LineWidth = 1.5;
       
%        q(i-2). Color(4) = 0.5;
%        q(i-2).LineWidth = 1.5;
    end
    if i>3
       p(i-3). Color(4) = 0.25;
       p(i-3).LineWidth = 1;
       
%        q(i-3). Color(4) = 0.25;
%        q(i-3).LineWidth = 1;
    end
   p(i) = plot(x(:,i),smooth(y(:,i),20),'LineWidth',2,'Color',[i/nt,1 - i/nt,i/nt]);
%    q(i) = plot(x(:,i),smooth(y2(:,i),20),'LineWidth',2,'Color',[i/nt,1 - i/nt,i/nt]);
   set(gca,'fontsize', 14);
   xlabel(xl,'interpreter','latex','FontSize',34)
   ylabel(yl,'interpreter','latex','FontSize',34)
   axis(range)
%    set(gca, 'xdir', 'reverse')
    axis square
   set(gcf,'color','w');
   set(gca, 'yScale', 'log')
   set(gca, 'xScale', 'log')
   if i==1
       f=getframe(gcf);
       [im,cm] = rgb2ind(f.cdata,256);
       imwrite(im,cm,fn,'DelayTime',dt,'LoopCount',inf);
    else
        f=getframe(gcf);
        [im,cm] = rgb2ind(f.cdata,256);
        imwrite(im,cm,fn,'DelayTime',dt,'writemode','append');
    end
   
   pause(dt)
end
hold off

end

