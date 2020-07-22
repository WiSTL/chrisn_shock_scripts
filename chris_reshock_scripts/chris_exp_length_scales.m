function [d,ds,tau] = chris_exp_length_scales(need_dpxdx,eta_range)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   1 - mixing thickness
%   2 - displacment thickness
%   3 - Momentum Thickness
%   4 - energy thickness
%   5 - stress thickness, ds - h_dot
%   6 - dissipation thickness
%   7-9 - taylor scale
%
%   d = integral <> - <>^2 dz
%   ds = integral < - ^2> dz - d
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

deltan=[];

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
    
%     v = v/deltadot;
    rho = rho/rho0;
    
    vm = mean(v');
    
    delta = delta/dpxdx;
    
    rm2 = interp1(eta,rm/rho0,eta_m,'pchip',0);
    vm2 = interp1(eta,vm,eta_m,'pchip',0);
    dvm2 = chris_derivative(vm2,abs(eta_m(2)-eta_m(1)));
    
    d12 = interp1(eta,mean((C.*(1-C))'),eta_m,'pchip',0);
    d22 = interp1(eta,mean((rho.*v)'),eta_m,'pchip',0);
    d32 = interp1(eta,mean((rho.*v.^2)'),eta_m,'pchip',0);
    d42 = interp1(eta,mean((rho.*v.^3)'),eta_m,'pchip',0);

    d(1,i) = delta;
    d(2,i) = delta*trapz(eta_m,(rm2.*vm2));
    d(3,i) = delta*trapz(eta_m,(rm2.*vm2).*((vm2)));
    d(4,i) = delta*trapz(eta_m,(rm2.*vm2).*((vm2.^2)));
    d(5,i) = delta*(dvm2(eta_m==0)).^-1;
    d(6,i) = delta*(trapz(eta_m,dvm2.^2)).^-1;
    
    ds(1,i) = delta-delta*trapz(eta_m,d12);
    ds(2,i) = delta*trapz(eta_m,d22)-d(2,i);
    ds(3,i) = delta*trapz(eta_m,d32)-d(3,i);
    ds(4,i) = delta*trapz(eta_m,d42)-d(4,i);
    ds(5,i) = deltadot;
    
    t0 = h_0./hd_0;
    
    tau(i) = (1/1000)*t_rs/t0;
    
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

    [l_ux,l_uy,l_vx,l_vy,l_cx,l_cy,l_cvx,l_cvy] = chris_taylor_microscale(u_p,v_p,Cpm);
    luxn(i) = l_ux/dpxdx;
    luyn(i) = l_uy/dpxdx;
    lvxn(i) = l_vx/dpxdx;
    lvyn(i) = l_vy/dpxdx;
    lcxn(i) = l_cx/dpxdx;
    lcyn(i) = l_cy/dpxdx;
    lcvxn(i) = l_cvx/dpxdx;
    lcvyn(i) = l_cvy/dpxdx;
    

%     d(7,i) = (1/dpxdx)*(l_cx^2 + l_cy^2)^0.5;
    d(7,i) = ((l_cx^2)./(l_cy^2))^0.5;
    d(8,i) = (1/dpxdx)*(l_ux^2 + l_uy^2)^0.5;
    d(9,i) = (1/dpxdx)*(l_vx^2 + l_vy^2)^0.5;
    d(10,i) = lcvyn(i)/lcvxn(i);
    d(11,i) = lcyn(i)/lcxn(i);
    
    [d(12,i),~,d(13,i)] = chris_v(eta,v,C);
    
end

end

