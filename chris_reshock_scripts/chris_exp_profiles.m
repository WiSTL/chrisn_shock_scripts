function [up2n,upvpn,apupn,apvpn,vp2n,vp3n,upupvpn,upvpdvpn,vpvpdupn,cn,cpn,un,vn,xin,udcn,vdcn,t11n,t12n,t22n,bn,t11Rn,t12Rn,t22Rn,rvn,run,eresn,esgsn,b11Rn,b12Rn,b22Rn,q2n,eta_m,Re] = chris_exp_profiles(need_dpxdx,eta_range)

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
    
    area_a = medfilt2(v);
    area_b = medfilt2(v);
    area_a(etam>0)=nan;
    area_b(etam<0)=nan;
    deltadot = delta_dot;%nanmean(max(area_a'))-nanmean(min(area_b'));
    
%     vm = mean2(v);
%     v = v-vm;
    
    [P,T,D,e,dd,up2,up3,vp2,vp3,cp2,cp3,upvp] = fluc_equation_terms(u,v,C,rho,deltadot,delta/dpxdx,mu./rho,eta);
  
    [tau,qc,S,sigma,S_favre,sigma_favre,epsilon_v,T,D,u_favre,v_favre,up,l,Re_t,K_sgs,K_res,e_res,e_sgs] = LES_filtered(u,v,C,rho,mu,dpxdx,dx_les);
    [t11,t12,t22,u_f,v_f,K,e,Ret,Ci,csaa,csab,csbb,ce,q2] = LES_profile(u,v,C,rho,mu,dpxdx);
    [t11R,t12R,t22R,KR] = reynolds_profile(u,v);
    [up21,upvp1,vp21,vp31,upupvp,upvpdvp,vpvpdup,um,vm,dum,dvm,dupvp,apup,apvp] = reynolds_mean_profile(u/deltadot,v/deltadot,eta,rho);
    [b] = chris_density_self_correlation(rho);
    [Xi,dcdx,dcdy] = chris_scalar_dissipation(C,dpxdx);
    [Ru,Rv] = chris_NS_Re(u,v,rho,mu);
    
    tkkR = t11R+t22R;
    bii = (t11R./tkkR)-(1/3);
    bij = (t12R./tkkR)-(1/3);
    bjj = (t22R./tkkR)-(1/3);
        
    udc = u.*dcdy;
    vdc = v.*dcdx;
    
    C2 = interp1(eta,mean(C'),eta_m,'pchip',NaN);
    Cp2 = interp1(eta,mean((c_p.^2)'),eta_m,'pchip',NaN);
    u2 = interp1(eta,mean((u/deltadot)'),eta_m,'pchip',NaN);
    v2 = interp1(eta,mean((v/deltadot)'),eta_m,'pchip',NaN);
    
    up2 = interp1(eta,up2,eta_m,'pchip',NaN);
    vp2 = interp1(eta,vp2,eta_m,'pchip',NaN);
    upvp2 = interp1(eta,upvp,eta_m,'pchip',NaN);
    apvp2 = interp1(eta,apvp,eta_m,'pchip',NaN);
    apup2 = interp1(eta,apup,eta_m,'pchip',NaN);
    vp3 = interp1(eta,vp3,eta_m,'pchip',NaN);
    upupvp = interp1(eta,upupvp,eta_m,'pchip',NaN);
    upvpdvp = interp1(eta,upvpdvp,eta_m,'pchip',NaN);
    vpvpdup = interp1(eta,vpvpdup,eta_m,'pchip',NaN);
%     apd2up = interp1(eta,apd2up,eta_m,'pchip',NaN);
%     apd2vp = interp1(eta,apd2vp,eta_m,'pchip',NaN);
    
%     rpvp2 = interp1(eta,rpvp,eta_m,'pchip',NaN);
    
    xi2 = interp1(eta,mean(Xi'),eta_m,'pchip',NaN);
    udc2 = interp1(eta,mean(udc'),eta_m,'pchip',NaN);
    vdc2 = interp1(eta,mean(vdc'),eta_m,'pchip',NaN);
    t112 = interp1(eta,mean(tau(:,:,1)'),eta_m,'pchip',NaN);
    t122 = interp1(eta,mean(tau(:,:,2)'),eta_m,'pchip',NaN);
    t222 = interp1(eta,mean(tau(:,:,3)'),eta_m,'pchip',NaN);
    b2 = interp1(eta,b,eta_m,'pchip',NaN);
    t11R2 = interp1(eta,t11R,eta_m,'pchip',NaN);
    t12R2 = interp1(eta,t12R,eta_m,'pchip',NaN);
    t22R2 = interp1(eta,t22R,eta_m,'pchip',NaN);
    ru2 = interp1(eta,Ru,eta_m,'pchip',NaN);
    rv2 = interp1(eta,Rv,eta_m,'pchip',NaN);
    eres2 = interp1(eta,e_res,eta_m,'pchip',NaN);
    esgs2 = interp1(eta,e_sgs,eta_m,'pchip',NaN);
    b11R2 = interp1(eta,bii,eta_m,'pchip',NaN);
    b12R2 = interp1(eta,bij,eta_m,'pchip',NaN);
    b22R2 = interp1(eta,bjj,eta_m,'pchip',NaN);
    q22 = interp1(eta,q2,eta_m,'pchip',NaN);
    
    Re(i) = (delta/dpxdx)*deltadot./(mean2(mu./rho));
    
    cn(i,:) = C2;
    cpn(i,:) = Cp2;
    un(i,:) = u2;
    vn(i,:) = v2;
    
    up2n(i,:) = up2;
    vp2n(i,:) = vp2;
    upvpn(i,:) = upvp2;
    apvpn(i,:) = apvp2;
    apupn(i,:) = apup2;
    vp3n(i,:) = vp3;
    upupvpn(i,:) = upupvp;
    upvpdvpn(i,:) = upvpdvp;
    vpvpdupn(i,:) = vpvpdup;
    
%     apd2upn(i,:) = apd2up;
%     apd2vpn(i,:) = apd2vp;
    
%     rpvpn(i,:) = rpvp2;
    
    xin(i,:) = xi2;
    udcn(i,:) = udc2;
    vdcn(i,:) = vdc2;
    t11n(i,:) = t112;
    t12n(i,:) = t122;
    t22n(i,:) = t222;
    bn(i,:) = b2;
    t11Rn(i,:) = t11R2;
    t12Rn(i,:) = t12R2;
    t22Rn(i,:) = t22R2;
    rvn(i,:) = rv2;
    run(i,:) = ru2;
    eresn(i,:) = eres2;
    esgsn(i,:) = esgs2;
    b11Rn(i,:) = b11R2;
    b12Rn(i,:) = b12R2;
    b22Rn(i,:) = b22R2;
    q2n(i,:) = q22;
    
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

