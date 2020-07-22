function [Pn,Tn,Dn,en,ddn,termn,eta_m] = chris_exp_equation_terms(need_dpxdx,eta_range)

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

%     u = medfilt2(medfilt2(u,[3 3]));
%     v = medfilt2(medfilt2(v,[3 3]));
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

    [P,T,D,e,dd,term] = fluc_equation_terms(u,v,C,rho,deltadot,delta/dpxdx,mu./rho,eta);
        
    P12 = interp1(eta,P(:,1),eta_m,'pchip',NaN);
    T12 = interp1(eta,T(:,1),eta_m,'pchip',NaN);
    D12 = interp1(eta,D(:,1),eta_m,'pchip',NaN);
    e12 = interp1(eta,e(:,1),eta_m,'pchip',NaN);
    
    P22 = interp1(eta,P(:,2),eta_m,'pchip',NaN);
    T22 = interp1(eta,T(:,2),eta_m,'pchip',NaN);
    D22 = interp1(eta,D(:,2),eta_m,'pchip',NaN);
    e22 = interp1(eta,e(:,2),eta_m,'pchip',NaN);
    
    P32 = interp1(eta,P(:,3),eta_m,'pchip',NaN);
    T32 = interp1(eta,T(:,3),eta_m,'pchip',NaN);
    D32 = interp1(eta,D(:,3),eta_m,'pchip',NaN);
    e32 = interp1(eta,e(:,3),eta_m,'pchip',NaN);
    
    P42 = interp1(eta,P(:,4),eta_m,'pchip',NaN);
    T42 = interp1(eta,T(:,4),eta_m,'pchip',NaN);
    D42 = interp1(eta,D(:,4),eta_m,'pchip',NaN);
    e42 = interp1(eta,e(:,4),eta_m,'pchip',NaN);
    
    P52 = interp1(eta,P(:,5),eta_m,'pchip',NaN);
    T52 = interp1(eta,T(:,5),eta_m,'pchip',NaN);
    D52 = interp1(eta,D(:,5),eta_m,'pchip',NaN);
    e52 = interp1(eta,e(:,5),eta_m,'pchip',NaN);

    dd2 = interp1(eta,dd,eta_m,'pchip',NaN);
    term2 = interp1(eta,term,eta_m,'pchip',NaN);
    
    Pn(1,i,:) = P12;
    Tn(1,i,:) = T12;
    Dn(1,i,:) = D12;
    en(1,i,:) = e12;
    
    Pn(2,i,:) = P22;
    Tn(2,i,:) = T22;
    Dn(2,i,:) = D22;
    en(2,i,:) = e22;
    
    Pn(3,i,:) = P32;
    Tn(3,i,:) = T32;
    Dn(3,i,:) = D32;
    en(3,i,:) = e32;
    
    Pn(4,i,:) = P42;
    Tn(4,i,:) = T42;
    Dn(4,i,:) = D42;
    en(4,i,:) = e42;
    
    Pn(5,i,:) = P52;
    Tn(5,i,:) = T52;
    Dn(5,i,:) = D52;
    en(5,i,:) = e52;
    
    ddn(i,:) = dd2;
    termn(i,:) = term2;
    
end
%     figure(1)
%     chris_plot_mean_std(eta_m,cn');
%     title('mean Concentration Profile')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('$\xi$')
%     
%     figure(2);
%     chris_plot_mean_std(eta_m,cpn');
%     title('mean Concentration fluctuation Profile')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('|$\xi$''|')
%     
%     figure(3);
%     chris_plot_mean_std(eta_m,un');
%     title('mean u')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('u')
% 
%     figure(5);
%     chris_plot_mean_std(eta_m,eresn');
%     title('mean resolved turbulent dissipation')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('$\epsilon_{res}$')
%     
%     figure(6);
%     chris_plot_mean_std(eta_m,esgsn');
%     title('mean subgrid turbulent dissipation')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('$\epsilon_{sgs}$')
%     
%     figure(7);
%     chris_plot_mean_std(eta_m,xin');
%     title('mean dissipation')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('u''')
%     
%     figure(8);
%     chris_plot_mean_std(eta_m,udcn');
%     title('mean udc')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('udc')
%     
%     figure(9);
%     chris_plot_mean_std(eta_m,vdcn');
%     title('mean vdc')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('vdc')
%     
%     figure(10);
%     chris_plot_mean_std(eta_m,bn');
%     title('density self correlation')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('b')
%     
%     figure(11);
%     chris_plot_mean_std(eta_m,t11Rn');
%     title('reynolds profile stress')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('$t11 R [m^2s^{-2}]$')
%     
%     figure(12);
%     chris_plot_mean_std(eta_m,t12Rn');
%     title('reynolds profile stress')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('t12 R')
%     
%     figure(13);
%     chris_plot_mean_std(eta_m,t22Rn');
%     title('reynolds profile stress')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('t22 R')
%     
%     figure(14);
%     chris_plot_mean_std(eta_m,run');
%     title('Ru')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('Ru')
%     
%     figure(15);
%     chris_plot_mean_std(eta_m,rvn');
%     title('Rv')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('Rv')
end

