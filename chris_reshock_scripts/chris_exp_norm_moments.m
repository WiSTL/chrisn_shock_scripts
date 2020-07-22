function [L,M,N,Cm,Cp2,Um,Up2,Wm,Wp2,Lambda2] = chris_exp_norm_moments(need_dpxdx,eta_range)

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

eta_m = eta_range;
C_m = 0:0.01:1;

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
    [xm,etam] = meshgrid(x,eta);
    
    [delta_dot,v,vm] = chris_v(eta,v,C);
     
    area_a = medfilt2(v);
    area_b = medfilt2(v);
    area_a(etam>0)=nan;
    area_b(etam<0)=nan;
    deltadot = nanmean(max(area_a'))-nanmean(min(area_b'));
     
    u = u/deltadot;
    v = v/deltadot;
    
    c_p = C;
    u_p = u;
    w_p = v;
    
    cm = nanmean(C,2);
    um = nanmean(u,2);
    wm = nanmean(v,2);
    
    for j = 1:nx
        c_p(:,j) = C(:,j) - cm;
        u_p(:,j) = u(:,j) - um;
        w_p(:,j) = v(:,j) - wm;
    end
    
    cp2 = nanmean(c_p.^2,2);
    up2 = nanmean(u_p.^2,2);
    wp2 = nanmean(w_p.^2,2);
    
    dcz = medfilt2(chris_gradient(c_p,1/dpxdx,1/dpxdx));
    dwz = medfilt2(chris_gradient(w_p,1/dpxdx,1/dpxdx));
    
    for n = 0:6
        for m = 0:6
            L(n+1,m+1,i,:) = interp1(eta,nanmean((u_p.^n).*(w_p.^m),2)./sqrt((up2.^n).*(wp2.^m)),eta_m,'pchip',NaN);
            M(n+1,m+1,i,:) = interp1(eta,nanmean((c_p.^n).*(w_p.^m),2)./sqrt((cp2.^n).*(wp2.^m)),eta_m,'pchip',NaN);
            N(n+1,m+1,i,:) = interp1(eta,nanmean((c_p.^n).*(u_p.^m),2)./sqrt((cp2.^n).*(up2.^m)),eta_m,'pchip',NaN);
            Lambda2(n+1,m+1,i,:) = interp1(eta,nanmean((c_p.^n).*(w_p.^m),2)./nanmean((dcz.^n).*(dwz.^m),2),eta_m,'pchip',NaN);            
        end
    end    
    
    Cm(i,:) = interp1(eta,cm,eta_m,'pchip',NaN);
    Cp2(i,:) = interp1(eta,cp2,eta_m,'pchip',NaN);
    
    Um(i,:) = interp1(eta,um,eta_m,'pchip',NaN);
    Up2(i,:) = interp1(eta,up2,eta_m,'pchip',NaN);
    
    Wm(i,:) = interp1(eta,wm,eta_m,'pchip',NaN);
    Wp2(i,:) = interp1(eta,wp2,eta_m,'pchip',NaN);
    
end

end

