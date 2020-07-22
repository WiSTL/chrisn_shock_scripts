function [T_hat,E_hat,v_hat,u_hat,C_hat,uv_hat,vv_hat,kdm,eta_m] = chris_exp_1D_homo_fourier(need_dpxdx,eta_range)

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
kdm = 0:0.0005:1;

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
    
    [delta_dot,v,vm] = chris_v(eta,v,C);

    Cmm = 0.5*erfc(etam);
    Cpm = C-Cmm;
    cm = mean(C');
    
    c_p = C;
    
    for j = 1:nx
        c_p(:,j) = C(:,j) - cm';
    end
    
    area_a = medfilt2(v);
    area_b = medfilt2(v);
    area_a(etam>0)=nan;
    area_b(etam<0)=nan;
    deltadot = delta_dot;%nanmean(max(area_a'))-nanmean(min(area_b'));
  
    
    [T_h,E_h,v_h,u_h,C_h,uv_h,vv_h,kd] = chris_fourier_homogenous_dir(C,u,v,delta,eta);

    for j = 1:nx
    T_hat1(:,j) = interp1(eta,T_h(:,j),eta_m,'pchip',0);
    E_hat1(:,j) = interp1(eta,E_h(:,j),eta_m,'pchip',0);
    v_hat1(:,j) = interp1(eta,v_h(:,j),eta_m,'pchip',0);
    u_hat1(:,j) = interp1(eta,u_h(:,j),eta_m,'pchip',0);
    C_hat1(:,j) = interp1(eta,C_h(:,j),eta_m,'pchip',0);
    uv_hat1(:,j) = interp1(eta,uv_h(:,j),eta_m,'pchip',0);
    vv_hat1(:,j) = interp1(eta,vv_h(:,j),eta_m,'pchip',0);
    end
    
    ne = length(eta_m);
    
    for j = 1:ne
    T_hat(j,:,i) = interp1(kd,T_hat1(j,:),kdm,'pchip',NaN);
    E_hat(j,:,i) = interp1(kd,E_hat1(j,:),kdm,'pchip',NaN);
    v_hat(j,:,i) = interp1(kd,v_hat1(j,:),kdm,'pchip',NaN);
    u_hat(j,:,i) = interp1(kd,u_hat1(j,:),kdm,'pchip',NaN);
    C_hat(j,:,i) = interp1(kd,C_hat1(j,:),kdm,'pchip',NaN);
    uv_hat(j,:,i) = interp1(kd,uv_hat1(j,:),kdm,'pchip',NaN);
    vv_hat(j,:,i) = interp1(kd,vv_hat1(j,:),kdm,'pchip',NaN);
    end

    vars = {'T_hat1','E_hat1','v_hat1','u_hat1','C_hat1','uv_hat1','vv_hat1'};
    clear(vars{:})
    
end

T_hat = mean(T_hat,3);
E_hat = mean(E_hat,3);
v_hat = mean(v_hat,3);
u_hat = mean(u_hat,3);
C_hat = mean(C_hat,3);
uv_hat = mean(uv_hat,3);
vv_hat = mean(vv_hat,3);

end

