function [] = chris_exp_1D_spectra_mean()

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

k_m = [100:0.05:1.4*10^4];
 
for i = 1:fileno
    i
    
     if fileno>1
        file = fileLocations{i};
    else
        file = fileLocations;
     end
   
     load(file);
     
     [ny,nx] = size(C);
     
     C(C>1)=1;
     C(C<0)=0;
     
     vm = mean2(v);
     v = v-vm;
     
     Cm = mean(C');
    y = [1:ny];
    x = [1:nx];

    ft = fittype('0.5*erfc((4/sqrt(2*pi))*(x-b)/a)');
    f = fit(y(:),Cm(:), ft,'StartPoint',[200,90]);

    delta = f.a;
    y0 = f.b;

    eta = (4/sqrt(2*pi))*(y-y0)/delta;

    [xm,etam] = meshgrid(x,eta);

    Cmm = 0.5*erfc(etam);
    Cpm = C-Cmm;
    
    if i==15 || i==16 || i==17
        dpxdx = 2274.375*2;
    else
        dpxdx = 2274.375;
    end

    [rho,mu] = chris_rho(C);

    r1 = 4.6081;
    r2 = 7.622;
    
    rho = r2+(r1-r2)*C;
    
    E = rho.*(u.^2 + v.^2);
    
%     [~,dru2dx] = chris_gradient(rho.*u.*u,1/dpxdx,1/dpxdx);
%     [druvdy,druvdx] = chris_gradient(rho.*u.*v,1/dpxdx,1/dpxdx);
%     [drv2dy,~] = chris_gradient(rho.*v.*v,1/dpxdx,1/dpxdx);
%     [dr2u] = chris_laplacian_2D(u,1/dpxdx,1/dpxdx);
%     
%     [dudy,dudx] = chris_gradient(u,1/dpxdx,1/dpxdx);
%     [dvdy,dvdx] = chris_gradient(v,1/dpxdx,1/dpxdx);
%     
%     E = -mu.*(-(2/3)*(dudx+dvdy).^2 + 2*dudx.^2 + 2*dvdy.^2 +(dudy+dvdx).^2);
    
    [eey,key] = power_spectra_1D(E,1/dpxdx,0,1);
    [eex,kex] = power_spectra_1D(E',1/dpxdx,0,1);
    
    [ecy,kcy] = power_spectra_1D(C,1/dpxdx,1,1);
    [ecx,kcx] = power_spectra_1D(C',1/dpxdx,1,1);
    
    aey = mean(eey)/max(mean(eey));
    aex = mean(eex)/max(mean(eex));
    acy = mean(ecy)/max(mean(ecy));
    acx = mean(ecx)/max(mean(ecx));
    
    aex2 = interp1(kex,aex,k_m,'pchip',NaN);
    aey2 = interp1(key,aey,k_m,'pchip',NaN);
    acx2 = interp1(kcx,acx,k_m,'pchip',NaN);
    acy2 = interp1(kcy,acy,k_m,'pchip',NaN);
    
    aexn(i,:) = aex2;
    acxn(i,:) = acx2;
    aeyn(i,:) = aey2;
    acyn(i,:) = acy2;
end

    figure(1);
    a = nanmean(aeyn);
    loglog(k_m,smooth(a/max(a)),'k','Linewidth',2)
    title('Energy Spectra in y direction')
    grid on
    hold on
    xlabel('k_y')
    ylabel('E_k(k_y)')
    
    figure(2)
    a = nanmean(aexn);
    loglog(k_m,smooth(a/max(a)),'k','Linewidth',2)
    title('Energy Spectra in x direction')
    grid on
    hold on
    xlabel('k_x')
    ylabel('E_k(k_x)')
    
    figure(3);
    a = nanmean(acyn);
    loglog(k_m,smooth(a/max(a)),'k','Linewidth',2)
    title('Scalar power Spectra in y direction')
    grid on
    hold on
    xlabel('k_y')
    ylabel('E_{\xi}(k_y)')

    figure(4)
    a = nanmean(acxn);
    loglog(k_m,smooth(a/max(a)),'k','Linewidth',2)
    grid on
    hold on
    title('Scalar power Spectra in x direction')
    xlabel('k_x')
    ylabel('E_{\xi}(k_x)')



end

