function [zeta_c_n,zeta_u_n,zeta_w_n,Sc2,Su2,Sw2,rx,Scn,Sun,Swn] = chris_exp_str_fn(need_dpxdx)

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
    
    [Sc2] = chris_no_structure_fn(C,2,dpxdx);
    [Su2] = chris_no_structure_fn(u,2,dpxdx);
    [Sw2,rx] = chris_no_structure_fn(v,2,dpxdx);
    
    k = log(squeeze(rx(1,2:4)))';
    
    Scn=[];
    Sun=[];
    Swn=[];
    
    for n=2:10
            [Scn(:,:,n)] = chris_no_structure_fn(C,n,dpxdx);
            [Sun(:,:,n)] = chris_no_structure_fn(u,n,dpxdx);
            [Swn(:,:,n)] = chris_no_structure_fn(v,n,dpxdx);
            
            a = log(squeeze(nanmean(Scn(:,2:4,n),1)))';
            a(isnan(a)) = 0;
            a(isinf(a)) = 0;
            f = fit(k,a, ft,'StartPoint',[1,3]);
            zeta_c_n(n,i) = f.b;
            
            a = log(squeeze(nanmean(Sun(:,2:4,n),1)))';
            a(isnan(a)) = 0;
            a(isinf(a)) = 0;
            f = fit(k,a, ft,'StartPoint',[1,3]);
            zeta_u_n(n,i) = f.b;
            
            a = log(squeeze(nanmean(Swn(:,2:4,n),1)))';
            a(isnan(a)) = 0;
            a(isinf(a)) = 0;
            f = fit(k,a, ft,'StartPoint',[1,3]);
            zeta_w_n(n,i) = f.b;
    end
    
    
%     z = [1:ny]/dpxdx;
%     
%     [K,E] = meshgrid(k,eta);
%     
%     E_c = E_c(:,1:nx/2);
%     E_c_e(:,:,i) = interp2(K,E,E_c,Kki,Ei);
%     
%     E_u = E_u(:,1:nx/2);
%     E_u_e(:,:,i) = interp2(K,E,E_u,Kki,Ei);
%     
%     E_w = E_w(:,1:nx/2);
%     E_w_e(:,:,i) = interp2(K,E,E_w,Kki,Ei);
    
    
    
end

end

