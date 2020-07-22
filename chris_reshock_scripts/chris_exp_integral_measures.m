function [deltan,hn,hmn,eps,scalar_eps,At,Re2,Retaylorc,Retayloru,Retaylorv,tkk_max,tkk,l_tc,l_tu,l_tv,ch,mh,tkk1,dav,hd0,h0,trs,ts] = chris_exp_integral_measures(need_dpxdx,eta_range)

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
hn=[];
hmn=[];
ch=[];
deltadotn=[];
deltadot2n=[];
Retn=[];
tkk=[];
Retaylor=[];
varc=[];
l_n=[];
l_t=[];
upn=[];
kn=[];
vmn=[];
luxn=[];
luyn=[];
lvxn=[];
lvyn=[];
lcxn=[];
lcyn=[];

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
    
    Cmm = 0.5*erfc(-etam);
    Cpm = C-Cmm;
    
    deltadot = 1;%delta_dot;

     [P,T,D,e,dd,up2,up3,vp2,vp3,cp2,cp3,apvp] = fluc_equation_terms(u,v,C,rho,deltadot,delta,mu./rho,eta);
    
    
    [dcdy,dcdx] = chris_gradient(C,1,1);
    
%     [C2] = shift_C(C);
%     
%     Cm = mean(C2');
    
    Cm = mean(C');
    
    vm = mean2(v);
    v = v+vm;
    
    vmn(i) = vm;

    b = C.*(1-C);

    deltadotn(i) = trapz(mean((v.*dcdy)').*(2*Cm-1)+mean((u.*dcdx)').*(2*Cm-1));

    deltadot2n(i) = trapz(4*mean(((0.5-C).*b.*(v-mean2(v)))'));
    
    [drvc] = chris_derivative(mean((rho.*v.*C)'),abs(eta(2)-eta(1)));
    [drv] = chris_derivative(mean((rho.*v)'),abs(eta(2)-eta(1)));
    rc = mean((rho.*C)');
    r = mean(rho');
    deltadot3n(i) = trapz(eta,(1./r).*(drvc+(rc./r).*drv).*(2*(rc./r)-1));
    
    area_a = v;
    area_b = v;
    area_a(etam<0)=nan;
    area_b(etam>0)=nan;
    deltadot4n(i) = delta_dot;%nanmean(max(area_a'))-nanmean(min(area_b'));
    
    deltan(i) = delta/dpxdx;
    
    [tau,qc,S,sigma,S_favre,sigma_favre,epsilon_v,T,D,u_favre,v_favre,up,l,Re_t,K_sgs,K_res,e_res,e_sgs,rho_gauss] = LES_filtered(u,v,C,rho,mu,dpxdx,10);
   
    a = chris_derivative(smooth(mean((rho.*(v-mean2(v)).*C)')),abs(eta(2)-eta(1)));
    
    deltadot3n(i) = 2*trapz(eta,mean((rho.*C)').*(a')./(mean(rho')).^2);
    
    Retn(i) = Re_t;
%     Re(i) = abs(delta*deltadotn(i)*mean2(rho)/mean2(mu));
    
%     Cpm(Cpm<=0)=NaN;
    
    varc(i) = nanmean(nanmean(abs(Cpm')));
    l_n(i) = l;
    upn(i)=mean(up);
    
    [t11,t12,t22,u_f,v_f,K,e,Ret] = LES_profile(u,v,C,rho,mu,dpxdx);
    [t11,t12,t22,up,vp,K] = reynolds_profile(u,v);
    
    
%   [h,i5,i95] = concentration_centerline_5_95(C);
%     hn(i) = abs(i5-i95);

    v_p = v;
    u_p = u;
    
    vm = mean(v');
    um = mean(u');
    
    for j = 1:nx
        v_p(:,j) = v(:,j) - vm';
    end
    
    for j = 1:nx
        u_p(:,j) = u(:,j) - um';
    end

    [l_ux,l_uy,l_vx,l_vy,l_cx,l_cy] = chris_taylor_microscale(u_p,v_p,Cpm);
    luxn(i) = l_ux/dpxdx;
    luyn(i) = l_uy/dpxdx;
    lvxn(i) = l_vx/dpxdx;
    lvyn(i) = l_vy/dpxdx;
    lcxn(i) = l_cx/dpxdx;
    lcyn(i) = l_cy/dpxdx;
%     
    vrms = sqrt(mean((vp2)));
    urms = sqrt(mean((up2)));
    l_tc(i) = (1/dpxdx)*(l_cx^2 + l_cy^2)^0.5;
    l_tu(i) = (1/dpxdx)*(l_ux^2 + l_uy^2)^0.5;
    l_tv(i) = (1/dpxdx)*(l_vx^2 + l_vy^2)^0.5;
%     
    Retaylorc(i) = mean(sqrt(vrms.^2+urms.^2))*l_tc(i)/nu_m;
    Retayloru(i) = mean(urms)*l_tu(i)/nu_m;
    Retaylorv(i) = mean(vrms)*l_tv(i)/nu_m;
%     
%     Re(i) = (3/20)*Retaylor(i)^2;
      Re2(i) = (delta/dpxdx)*delta_dot./nu_m;
      ch(i) = (delta_dot*dpxdx/delta).*(h_0/hd_0);
      mh(i) = delta_dot/1000;
%     
%     [~,dru2] = chris_gradient(rho.*u.*u,1,1);
%     [drv2,~] = chris_gradient(rho.*v.*v,1,1);
%     [druvdy,druvdx] = chris_gradient(rho.*u.*v,1,1);
% 
%     d2u = chris_laplacian_2D(u,1,1);
%     d2v = chris_laplacian_2D(v,1,1);
% 
%     Ru(i) = mean(sqrt(mean(((dru2+druvdy).^2)')./mean(((mu.*d2u).^2)')));
%     Rv(i) = mean(sqrt(mean(((drv2+druvdx).^2)')./mean(((mu.*d2v).^2)')));
%     
      dpx(i) = dpxdx;
%     
      tkki = up2+vp2;
%     
%     vm = mean(v');
%     
%     v2 = interp1(eta,vm,eta_m,'pchip',0);

      eta_tkk = [-0.5:0.01:0.5];
      tkk2 = interp1(eta,tkki,eta_tkk,'pchip',0);
%     P2 = interp1(eta,P(:,2),eta_m,'pchip',0);
%     T2 = interp1(eta,T(:,2),eta_m,'pchip',0);
%     eps2 = interp1(eta,eps,eta_m,'pchip',0);
      cp22 = interp1(eta,cp2,eta_tkk,'pchip',0);
      up22 = interp1(eta,up2,eta_tkk,'pchip',0);
      vp22 = interp1(eta,vp2,eta_tkk,'pchip',0);
      apvp2 = interp1(eta,apvp,eta_tkk,'pchip',0);
%     
%     
%     
      tkk(i) = trapz(eta_tkk,tkk2)./delta_dot.^2;
%     Pi(i) = trapz(eta_m,P2);
%     Ti(i) = trapz(eta_m,T2);
%     ei(i) = trapz(eta_m,eps2);

     [hmn(i),dv(i),dav(i)] = chris_fluc_length_scales(eta_tkk,cp22,up22,vp22,apvp2);
%       hn(i) = trapz(eta_m,cp22);
%       du(i) = trapz(eta_m,up22);
%       dv(i) = trapz(eta_m,vp22);
%     hmn(i) = trapz(eta_m,(0.5-v2).*(0.5+v2)/deltadot^2);
%     
%     tkk_units(i) = delta*tkk(i)*deltadot4n(i)^2;
      hn(i) = trapz(eta_tkk,cp22);
      tkk_max(i) = nanmax(tkk2);
      tkk1(i) = dv(i);

    trs(i) = t_rs;
    ts(i) = t_s;
    hd0(i) = hd_0;
    h0(i) = h_0;
    At(i) = A_rs;
    
    cmean = mean(C');
    c_p = C;
    
    
    for j = 1:nx
        c_p(:,j) = C(:,j) - cmean';
    end
    
    cp2 = mean(c_p(:).^2);
    
    [Svx,rx] = chris_no_structure_fn(v,3,dpxdx);
    eps(i) = nanmean(nanmean((5/4)*mean(Svx(:,2:4),1)./rx(2:4)))*(delta/dpxdx)/delta_dot^3;
    
    [Scx,rx] = chris_no_structure_fn(C,3,dpxdx);   
    scalar_eps(i) = nanmean(nanmean((5/4)*mean(Scx(:,2:4),1)./rx(2:4)))*(delta/dpxdx)/(delta_dot*cp2);
    
end

% figure(1)
% plot(deltan,hn,'*')
% 
% figure(2)
% plot(deltadot3n,deltadot4n,'*')
% 
% figure(3)
% plot(Retn,varc,'*')
% 
% figure(4)
% plot(Re,l_n,'*')
% 
% figure(5)
% plot(Retn,upn,'*')
% 
% figure(6)
% % plot(Retn,deltan,'+')
% % hold on
% plot(Re,deltan,'+')
% hold on
% plot(sqrt(Rv.^2+Ru.^2),deltan,'+')
% % plot(Ru,deltan,'+')
% % plot(Rv,deltan,'+')
% 
% figure(9)
% plot(deltan,vmn,'*')
% 
% figure(10)
% plot(deltadotn,deltadot3n,'*')
% 
% figure(11)
% plot(abs(luxn)*1000,abs(luyn)*1000,'*')
% hold on
% 
% figure(12)
% plot(abs(lvxn)*1000,abs(lvyn)*1000,'*')
% hold on
% 
% figure(13)
% plot(abs(lcxn)*1000,abs(lcyn)*1000,'*')
% hold on
% 
% figure(14)
% plot(Re,Retaylor,'+')
% 
% figure(15)
% plot(Re,Ru,'+')
% hold on
% plot(Re,Rv,'+')
% 
% figure(16)
% plot(Re,sqrt(Rv.^2+Ru.^2),'+')
% 
figure(17)
plot(abs(luxn)*1000,abs(luyn)*1000,'*')
hold on

figure(18)
plot(abs(lvxn)*1000,abs(lvyn)*1000,'*')
hold on

figure(19)
plot(abs(lcxn)*1000,abs(lcyn)*1000,'*')
hold on

figure(20)
plot(abs(lvxn)*1000,abs(luxn)*1000,'*')
hold on

figure(21)
plot(abs(lvyn)*1000,abs(luyn)*1000,'*')
hold on
% 
% figure(20)
% plot(10*(lcxn./(deltan.*dpx)).^(-2),10*(lcyn./(deltan.*dpx)).^(-2),'*');
% 
% figure(22)
% plot(deltan,tkk);
% 
% figure(23)
% plot(abs(lvxn./(deltan.*dpx)),abs(luxn./(deltan.*dpx)),'*')
% 
% figure(24)
% plot(abs(lvyn./(deltan.*dpx)),abs(luyn./(deltan.*dpx)),'*')

% figure(25)
% plot(trs,deltan,'*')
% 
% figure(26)
% plot(trs,deltadot4n./hd0,'*')

end

