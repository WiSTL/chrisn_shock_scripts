function [bs,cs,ds,es] = chris_exp_mean_profiles(need_dpxdx)

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
     
     [delta_dot,v,vm] = chris_v(eta,v,C);

    [xm,etam] = meshgrid(x,eta);

    Cmm = 0.5*erfc(etam);
    Cpm = C-Cmm;
    
    Cm = mean(C');
    rm = mean(rho');
    
    C_p = C;
    v_p = v;
    r_p = rho;
    
    for j = 1:nx
        C_p(:,j) = C(:,j) - Cm';
    end
    
    for j = 1:nx
        r_p(:,j) = rho(:,j) - rm';
    end
    
    drm = chris_derivative(smooth(rm,20),abs(eta(1)-eta(2)));
    d2rm = chris_derivative(smooth(drm,20),abs(eta(1)-eta(2)));
    
    rp2 = mean((r_p.^2)');
    rp3 = mean((r_p.^3)');
    rp4 = mean((r_p.^4)');
    drp2 = chris_derivative(smooth(rp2,20),abs(eta(1)-eta(2)));
    
    
    area_a = medfilt2(v);
    area_b = medfilt2(v);
    area_a(etam>0)=nan;
    area_b(etam<0)=nan;
    deltadot = abs(2*nanmean(min(area_b')));

    [P,T,D,e,dd,up2,vp2,upvp] = fluc_equation_terms(u,v,C,rho,delta_dot,delta/dpxdx,mu./rho,eta);
    
    [tau,qc,S,sigma,S_favre,sigma_favre,epsilon_v,T,D,u_favre,v_favre,up,l,Re_t,K_sgs,K_res,e_res,e_sgs] = LES_filtered(u,v,C,rho,mu,dpxdx,dx_les);
    [t11,t12,t22,u_f,v_f,K,e,Ret] = LES_profile(u,v,C,rho,mu,dpxdx);
    [t11R,t12R,t22R,KR] = reynolds_profile(u,v);
    [b] = chris_density_self_correlation(rho);
    [Xi,dcdx,dcdy] = chris_scalar_dissipation(C,dpxdx);
    [Ru,Rv] = chris_NS_Re(u,v,rho,mu);
    
    
    u = u/delta_dot;
    v = v/delta_dot;
    
    udc = u.*dcdy;
    vdc = v.*dcdx;
    
    vm2 = mean(v');
    
    for j = 1:nx
        v_p(:,j) = v(:,j) - vm2';
    end
    
    a = chris_derivative(smooth(mean((rho.*v.*C)'),10),abs(eta(2)-eta(1)));
    
    vmdcm = smooth((smooth(vm2,20)').*chris_derivative(smooth(Cm,40),abs(eta(2)-eta(1)))',20)';
    
    [dcp,~] = chris_gradient(medfilt2(C_p),1/dpxdx,1/dpxdx);
    
    cpdvp = smooth(mean((medfilt2(v_p.*dcp))'),20)';
    
    dcpvp = chris_derivative(smooth(mean((medfilt2(v_p.*C_p))'),20)',abs(eta(2)-eta(1)));
    
    b = (2*Cm-1).*(vmdcm);
    
    c = (2*Cm-1).*cpdvp;
    
    d = dcpvp;
    
    bs(i) = (delta/dpxdx)*trapz(eta,b);
    cs(i) = (delta/dpxdx)*trapz(eta,c);
    ds(i) = delta_dot*delta/mean2(mu./rho);
    es(i) = deltadot*delta/mean2(mu./rho);
    
%     figure(1);
%     plot(eta,mean(C'))
%     title('mean Concentration Profile')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('$\xi$')
%     
%     figure(2);
%     plot(eta,mean(abs((C_p.^2)')))
%     title('mean Concentration fluctuation Profile')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('|$\xi$''|')
%     
%     figure(3);
%     plot(eta,mean(u'))
%     title('mean u')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('u')
% 
%     figure(5);
%     plot(eta,mean(v'))
%     title('mean resolved turbulent dissipation')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('mean v')
%     
%     figure(6);
%     plot(eta,up2)
%     title('u variance')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('up2')
%     
%     figure(7);
%     plot(eta,vp2)
%     title('v variance')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('vp2')
%     
%     figure(8);
%     plot(eta,vmdcm)
%     title('rm.*drm')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('rm.*drm')
%     
%     figure(9);
%     plot(eta,cpdvp)
%     title('rm.*d2rm')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('rm.*d2rm')
%     
%     figure(10);
%     plot(eta,smooth(b))
%     title('rm.*d2rm')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('rm.*d2rm')
%     
%     figure(11);
%     plot(eta,smooth(c))
%     title('rm.*d2rm')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('rm.*d2rm')
%     
%     figure(12);
%     plot(eta,smooth(d))
%     title('rm.*d2rm')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('rm.*d2rm')
    
    
%     figure(6);
%     plot(eta,e_sgs)
%     title('mean subgrid turbulent dissipation')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('$\epsilon_{sgs}$')
%     
%     figure(7);
%     plot(eta,mean(Xi'))
%     title('mean dissipation')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('u''')
%     
%     figure(8);
%     plot(eta,mean(udc'))
%     title('mean udc')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('udc')
%     
%     figure(9);
%     plot(eta,mean(vdc'))
%     title('mean vdc')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('vdc')
%     
%     figure(10);
%     plot(eta,b)
%     title('density self correlation')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('b')
%     
%     figure(11);
%     plot(eta,t11R)
%     title('reynolds profile stress')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('t11 R [m^2s^{-2}]')
%     
%     figure(12);
%     plot(eta,t12R)
%     title('reynolds profile stress')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('t12 R')
%     
%     figure(13);
%     plot(eta,t22R)
%     title('reynolds profile stress')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('t22 R')
%     
%     figure(14);
%     plot(eta,smooth(abs(Ru)'));
%     title('Ru')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('Ru')
%     
%     figure(15);
%     plot(eta,smooth(abs(Rv)'));
%     title('Rv')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('Rv')
%     
%     figure(16);
%     plot(eta,a);
%     title('derivative of <\rho*v*C>')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('d<\rho*v*C>/d\eta')
%     
%     figure(17);
%     plot(eta,up2);
%     title('up2')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('up2')
%     
%     figure(18);
%     plot(eta,vp2);
%     title('vp2')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('vp2')
%     
%     figure(19);
%     plot(eta,upvp);
%     title('upvp')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('upvp')
%     
%     figure(20);
%     plot(eta,up2+vp2);
%     title('tkk')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('tkk')
%     
%     b = 4*mean(C').*(1-mean(C'));
%     
%     figure(21);
%     plot(eta,b)
%     title('b')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('$b$')
%     
%     figure(22);
%     plot(delta,trapz(eta,b))
%     title('b')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('$b$')
%     
%     figure(23);
%     plot(eta,cdvc)
%     title('b')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('$cdvc$')
%     
%     figure(24);
%     plot(eta,cdvpcp)
%     title('b')
%     grid on
%     hold on
%     xlabel('$\eta$')
%     ylabel('$cdvpcp$')
    
    
end

end

