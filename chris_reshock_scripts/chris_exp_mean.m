function [cnn,unn,vnn,e,x] = chris_exp_mean(need_dpxdx,eta_range,x_range)

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
x_m = x_range;

for i = 1:fileno
    i
    
    clear c_p;
    
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
    cm = mean(C');
    
    for j = 1:nx
        c_p(:,j) = C(:,j) - cm';
    end
    
    area_a = medfilt2(v);
    area_b = medfilt2(v);
    area_a(etam>0)=nan;
    area_b(etam<0)=nan;
    deltadot = nanmean(max(area_a'))-nanmean(min(area_b'));
    
    vm = mean2(v);
    v = (v-vm)/deltadot;
    
    um = mean2(u);
    u = (u-um)/deltadot;
    
    for j = 1:nx
       cn(:,j) =  interp1(eta,C(:,j),eta_m,'pchip',0);
       un(:,j) =  interp1(eta,u(:,j),eta_m,'pchip',0);
       vn(:,j) =  interp1(eta,v(:,j),eta_m,'pchip',0);
    end
    
    [ny,nx] = size(cn);
    
    for j = 1:ny
         cnn(i,j,:) =  interp1([1:nx]/dpxdx,cn(j,:),x_m,'pchip',NaN);
         unn(i,j,:) =  interp1([1:nx]/dpxdx,un(j,:),x_m,'pchip',NaN);
         vnn(i,j,:) =  interp1([1:nx]/dpxdx,vn(j,:),x_m,'pchip',NaN);
    end
    
    [x,e] = meshgrid(x_m,eta_m);
    
end
    
    