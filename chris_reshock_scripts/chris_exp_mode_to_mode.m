function [Tuub,Tccb,Tuum,Tccm,Tuut,Tcct] = chris_exp_mode_to_mode(need_dpxdx)

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

deltan=[];

etai = -2:0.05:2;
Ki = 0.25*[0:5:2*10^3];

[Kki,Ei] = meshgrid(Ki,etai);

Ei = flipud(Ei);

ft = fittype('a+b*x');

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
     
     [C,eta,x,delta] = chris_C(C);
    
    delta = delta/dpxdx;
    
    t0 = h_0./hd_0;
    
    tau(i) = (1/1000)*t_rs/t0;
    
    [xm,etam] = meshgrid(x,eta);
    [xi,Ei] = meshgrid(x,etai);
    
    area_a = imgaussfilt(medfilt2(v));
    area_b = imgaussfilt(medfilt2(v));
    area_a(etam>0)=nan;
    area_b(etam<0)=nan;
    
    area_a(etam<-0.5)=nan;
    area_b(etam>0.5)=nan;
    
    deltadot = nanmean(max(area_a'))-nanmean(min(area_b'));
    
    u=medfilt2(u(:,1:nx),[3 3])/deltadot;
    v=medfilt2(v(:,1:nx),[3 3])/deltadot;
    C=medfilt2(C(:,1:nx),[3 3]);
    
    C = interp2(xm,etam,C,xi,Ei);
    u = interp2(xm,etam,u,xi,Ei);
    v = interp2(xm,etam,v,xi,Ei);
    
    u(isnan(u))=0;
    v(isnan(v))=0;
    C(isnan(C))=0;
    
    up = u - nanmean(u,2);
    vp = v - nanmean(v,2);
    cp = C - nanmean(C,2);
        
    [Tccb(:,:,i),Tuub(:,:,i)] = chris_mode_to_mode_transfer(up(11:31,:),vp(11:31,:),cp(11:31,:));
    [Tccm(:,:,i),Tuum(:,:,i)] = chris_mode_to_mode_transfer(up(31:51,:),vp(31:51,:),cp(31:51,:));
    [Tcct(:,:,i),Tuut(:,:,i)] = chris_mode_to_mode_transfer(up(51:71,:),vp(51:71,:),cp(51:71,:));
    
end

end

