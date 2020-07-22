
function [d,tau] = chris_exp_Taylor_Scales(need_dpxdx,eta_range)

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
     
     [C,~,~,delta] = chris_C(C);
    
    delta = delta/dpxdx;
    
    t0 = h_0./hd_0;
    
    tau(i) = (1/1000)*t_rs/t0;
    
    v_p=[];
    u_p=[];
    C_p=[];
    
    v_p = v;
    u_p = u;
    C_p = C;
    
    vm = mean(v');
    um = mean(u');
    Cm = mean(C');
    
    for j = 1:nx
        v_p(:,j) = v(:,j) - vm';
        u_p(:,j) = u(:,j) - um';
        C_p(:,j) = C(:,j) - Cm';
    end

    [l_ux,l_uy,l_vx,l_vy,l_cx,l_cy,l_cvx,l_cvy,l_cux,l_cuy,l_K,l_C] = chris_taylor_microscale(u_p,v_p,C_p,dpxdx);
    
    d(1,i) = l_ux/delta;
    d(2,i) = l_uy/delta;
    d(3,i) = l_vx/delta;
    d(4,i) = l_vy/delta;
    d(5,i) = l_cx/delta;
    d(6,i) = l_cy/delta;
    d(7,i) = l_cvx/delta^2;
    d(8,i) = l_cvy/delta^2;
    d(9,i) = l_cux/delta^2;
    d(10,i) = l_cuy/delta^2;
    
    vp2 = mean(v_p.^2,2);
    up2 = mean(u_p.^2,2);
    
    d(11,i) = trapz(smooth(vp2))./trapz(smooth(up2));
    
    for j = 1:ny
        [Ew(j,:),Kk(j,:)] = pspectrum(v_p(j,:),[1:nx]/dpxdx,'leakage',0.85);
        [Eu(j,:),Kk(j,:)] = pspectrum(u_p(j,:),[1:nx]/dpxdx,'leakage',0.85);
         Ek(j,:) = Ew(j,:)+Eu(j,:);
        [Ec(j,:),Kc(j,:)] = pspectrum(C_p(j,:),[1:nx]/dpxdx,'leakage',0.85);
        
        Lk(j,i) = trapz(squeeze(Kk(j,2:end)),squeeze(Ek(j,2:end))./squeeze(Kk(j,2:end)))./trapz(squeeze(Kk(j,2:end)),squeeze(Ek(j,2:end)));
        Lc(j,i) = trapz(squeeze(Kc(j,2:end)),squeeze(Ec(j,2:end))./squeeze(Kc(j,2:end)))./trapz(squeeze(Kc(j,2:end)),squeeze(Ec(j,2:end)));
    end
    
    d(12,i) = trapz(1/dpxdx,Lk(:,i))/delta;
    d(13,i) = trapz(1/dpxdx,Lc(:,i))/delta;
    
    d(14,i) = l_K/delta;
    d(15,i) = l_C/delta;
    
end

end

