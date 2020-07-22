function [un,up2n,up3n,vn,vp2n,vp3n,cn,cp2n,cp3n,upvp2n,luxn,luyn,lvxn,lvyn,lcxn,lcyn,epsxn,epsyn,Rt,Re2,del,beta,tau,eta_m] = chris_exp_3_moments(need_dpxdx,eta_range)

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
    
    vmm = mean(v');
    
    area_a = imgaussfilt(medfilt2(v));
    area_b = imgaussfilt(medfilt2(v));
    area_a(etam>0)=nan;
    area_b(etam<0)=nan;
    
    area_a(etam<-0.5)=nan;
    area_b(etam>0.5)=nan;
    
    dv = hd_0/(0.28*A_rs);
    
    deltadot = nanmean(max(area_a'))-nanmean(min(area_b'));
    
%     vm = mean2(v);
%     v = v-vm;
    
    v_p = v;
    u_p = u;
    
    vm = mean2(v');
    um = mean2(u');
    
    for j = 1:nx
        v_p(:,j) = v(:,j) - vm';
    end
    
    for j = 1:nx
        u_p(:,j) = u(:,j) - um';
    end
    
    v_p = v_p/deltadot;
    u_p = u_p/deltadot;
    
    [dupdy,dupdx] = chris_gradient(u_p,1/dpxdx,1/dpxdx);
    [dvpdy,dvpdx] = chris_gradient(v_p,1/dpxdx,1/dpxdx);
    
    epsx = 2*nu_m*(dupdx.^2);
    epsx = nanmean(medfilt2(epsx,[10 10])');
    
    epsy = 2*nu_m*(dvpdy.^2);
    epsy = nanmean(medfilt2(epsy,[10 10])');

    [l_ux,l_uy,l_vx,l_vy,l_cx,l_cy] = chris_taylor_microscale_spatial(u_p,v_p,c_p);

    [P,T,D,e,dd,up2,up3,vp2,vp3,cp2,cp3,upvp] = fluc_equation_terms(u_p,v_p,C,rho,deltadot,delta/dpxdx,mu./rho,eta);
    
    C2 = interp1(eta,mean(rho'),eta_m,'pchip',NaN);
    Cp2 = interp1(eta,cp2,eta_m,'pchip',NaN);
    Cp3 = interp1(eta,cp3,eta_m,'pchip',NaN);
    u2 = interp1(eta,mean((u/deltadot)'),eta_m,'pchip',NaN);
    v2 = interp1(eta,mean((v/deltadot)'),eta_m,'pchip',NaN);
    
    upvp2 = interp1(eta,upvp,eta_m,'pchip',NaN);
    up2 = interp1(eta,up2,eta_m,'pchip',NaN);
    up3 = interp1(eta,up3,eta_m,'pchip',NaN);
    vp2 = interp1(eta,vp2,eta_m,'pchip',NaN);
    vp3 = interp1(eta,vp3,eta_m,'pchip',NaN);
    
    lux2 = interp1(eta,l_ux/dpxdx,eta_m,'pchip',NaN);
    luy2 = interp1(eta,l_uy/dpxdx,eta_m,'pchip',NaN);
    lvx2 = interp1(eta,l_vx/dpxdx,eta_m,'pchip',NaN);
    lvy2 = interp1(eta,l_vy/dpxdx,eta_m,'pchip',NaN);
    lcx2 = interp1(eta,l_cx/dpxdx,eta_m,'pchip',NaN);
    lcy2 = interp1(eta,l_cy/dpxdx,eta_m,'pchip',NaN);
    
    epsx2 = interp1(eta,epsx,eta_m,'pchip',NaN);
    epsy2 = interp1(eta,epsy,eta_m,'pchip',NaN);
    
    lu = sqrt(lux2.^2+luy2.^2);
    lv = sqrt(lvx2.^2+lvy2.^2);
    lambda = sqrt(lu.^2+lv.^2);
    
    Rt(i,:) = sqrt(up2+vp2).*lambda./nu_m;
    Re2(i) = (delta/dpxdx)*delta_dot./nu_m;
    
    Rt(i,:) = Rt(i,:)./sqrt(Re2(i));
    
    if exist('nu_m','var')
    Re(i) = (delta/dpxdx)*deltadot./nu_m;
    else
    Re(i) = (delta/dpxdx)*deltadot./(mean2(mu./rho)); 
    end
    t_0 = h_0./hd_0;
    tau(i) = (1/1000).*t_rs/t_0;
    
    cn(i,:) = C2;
    cp2n(i,:) = Cp2;
    cp3n(i,:) = Cp3;
    un(i,:) = u2;
    vn(i,:) = v2;
    
    upvp2n(i,:) = upvp2;
    up2n(i,:) = up2;
    up3n(i,:) = up3;
    vp2n(i,:) = vp2;
    vp3n(i,:) = vp3;
    
    luxn(i,:) = lux2;
    luyn(i,:) = luy2;
    lvxn(i,:) = lvx2;
    lvyn(i,:) = lvy2;
    lcxn(i,:) = lcx2;
    lcyn(i,:) = lcy2;
    
    epsxn(i,:) = epsx2;
    epsyn(i,:) = epsy2;
    
    del(i,:) = dpxdx*(lvy2 - lux2)/delta;
    beta(i,:) = vp2./up2;
    
    
    if exist('nu_m','var')
        clear nu_m
    end
    
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

