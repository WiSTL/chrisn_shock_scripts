function [un,vn] = chris_exp_contour_averages(need_dpxdx)

[filename,pathname] = uigetfile('*.mat','multiselect','on');

fileno = size(filename);
fileno=fileno(2)

if fileno>1
    for i = 1:fileno
        fileLocations{i} = [pathname filename{i}];
    end
    
else
    
    fileLocations = [pathname filename];
    
end

un=[];
vn=[];
 
for i = 1:fileno
    i
    
     if fileno>1
        file = fileLocations{i};
    else
        file = fileLocations;
     end
   
     load(file);
     
     dx_les = 10;
     if need_dpxdx
        if i==15 || i==16 || i==17
            dpxdx = 2274.375*2;
            dx_les = 10;
        else
            dpxdx = 2274.375;
            dx_les = 5;
        end
     end
     [ny,nx] = size(C);
     [C,eta,x,delta,y0] = chris_C(C);
     [rho,mu] = chris_rho(C);

    [xm,etam] = meshgrid(x,eta);

    Cmm = 0.5*erfc(etam);
    Cpm = C-Cmm;
    
    area_a = medfilt2(v);
    area_b = medfilt2(v);
    area_a(etam>0)=nan;
    area_b(etam<0)=nan;
    deltadot = nanmean(max(area_a'))-nanmean(min(area_b'));
    
    us = (u-nanmean(nanmean(u)))/deltadot;
    vs = (v-nanmean(nanmean(v)))/deltadot;
    
    n=0;
    for j = 0:0.001:1
        n=n+1;
        
        ch = C(C<=j+0.0005);
        uh = us(C<=j+0.0005);
        vh = vs(C<=j+0.0005);
        
        un(i,n) = nanmean(uh(ch>=j-0.0005));
        vn(i,n) = nanmean(vh(ch>=j-0.0005));
    end
    
    
end


end

