function [beta,d,ch,reh,ref,omega,tau] = chris_exp_anisotropy(need_dpxdx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   beta - w'2/u'2
%   delta - lambda_w - lambda_u
%   omega - Re_t/sqrt(Re_h)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

eta_m = -0.5:0.01:0.5;
 
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

     rho0 = 7.622;
     
    [xm,etam] = meshgrid(x,eta);

    Cmm = 0.5*erfc(etam);
    Cpm = C-Cmm;
    rm = mean(rho');
    
    area_a = medfilt2(v);
    area_b = medfilt2(v);
    area_a(etam>0)=nan;
    area_b(etam<0)=nan;
    deltadot = nanmean(max(area_a'))-nanmean(min(area_b'));
    
    [dcdy,dcdx] = chris_gradient(C,1,1);
    
    Cm = mean(C');
    
    vm = mean(v');
    
    t0 = h_0./hd_0;
    tau(i) = (1/1000)*t_rs/t0;
    
    v_p = v;
    u_p = u;
    
    vm = mean2(v');
    um = mean(u');
    
    for j = 1:nx
        v_p(:,j) = v(:,j) - vm';
    end
    
    for j = 1:nx
        u_p(:,j) = u(:,j) - um';
    end
    
    up2 = mean((u_p.^2)');
    vp2 = mean((v_p.^2)');
    
    vp2 = interp1(eta,vp2,eta_m,'pchip',0);
    up2 = interp1(eta,up2,eta_m,'pchip',0);
    
    up2 = delta*trapz(eta_m,up2);
    vp2 = delta*trapz(eta_m,vp2);
    
    delta = delta/dpxdx;
    
    [l_ux,l_uy,l_vx,l_vy,l_cx,l_cy] = chris_taylor_microscale(u_p,v_p,Cpm);
    luxn(i) = l_ux/dpxdx;
    luyn(i) = l_uy/dpxdx;
    lvxn(i) = l_vx/dpxdx;
    lvyn(i) = l_vy/dpxdx;
    lcxn(i) = l_cx/dpxdx;
    lcyn(i) = l_cy/dpxdx;
    
    ret = sqrt(up2+vp2)*lvyn(i)/(nu_m);
    
    beta(i) = vp2/up2;
    d(i) = (lvyn(i) - luxn(i))/delta;
    reh(i) = delta*deltadot/nu_m;
    omega(i) = ret/sqrt(reh(i));
    ch(i)  = (deltadot/delta)*(h_0./hd_0);
    
    dv = hd_0/(0.28*A_rs);
    
    ref(i) = dv*h_0/nu_m;
    
end

end

