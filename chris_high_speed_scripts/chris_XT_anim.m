function  [p] = chris_XT_anim(t,y0,h,range,dt,fn)

nt = length(t);

figure
hold on

for i = 1:nt
        
    y1 = y0(1:i);
    y2 = y0(1:i)+h(1:i)/2;
    y3 = y0(1:i)-h(1:i)/2;

    X=[t(1:i),fliplr(t(1:i))];                %#create continuous x value array for plotting
    Y=[y2,fliplr(y3)];              %#create y values for out and then back
    m = fill(X,Y,[i/nt,1 - i/nt,i/nt]); 

    set(m,'FaceAlpha',0.3);
    set(m,'edgealpha',0)

    hold on;
    p(i) = plot(t(1:i),y1,'color',[i/nt,1 - i/nt,i/nt],'linewidth',2);
    hold off

   set(gca,'fontsize', 14);
   xlabel('$t_s\ [\rm{s}]$','interpreter','latex','FontSize',34)
   ylabel('$\rm{Z} [\rm{m}]$','interpreter','latex','FontSize',34)
   axis(range)
   set(gca, 'ydir', 'reverse')
    axis square
   set(gcf,'color','w');
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

