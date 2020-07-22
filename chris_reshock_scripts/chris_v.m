function [delta_dot,v,vm,cv] = chris_v(eta,v,C)

    [ny,nx] = size(C);
    x = [1:nx];
    
    [xm,etam] = meshgrid(x,eta);

    vm = mean2(v);
    
    v = v-vm;
    
    area_a = imgaussfilt(medfilt2(v));
    area_b = imgaussfilt(medfilt2(v));
    area_a(etam<0)=nan;
    area_b(etam>0)=nan;
    
    area_a(etam>0.5)=nan;
    area_b(etam<-0.5)=nan;
    
    delta_dot = nanmean(max(area_a'))-nanmean(min(area_b'));

    cm = mean(C');
    
    for i = 1:nx
        Cp(:,i) = C(:,i) - cm';
    end
    
    itop = find(eta>-0.2);
    range = find(eta(itop)<0.2);
    range = itop(1):itop(range(end));
    
    cv = mean((Cp.*v)');
%     delta_dot = -4*mean(cv(range));
    
end

