function [E_c_e,Pi_w_e,Pi_u_e,Pi_c_e,Pi_z_w_e,Pi_z_u_e,Pi_z_c_e,P_c_e,P_u_e,P_w_e,E_u_e,E_w_e,D_c_e,D_u_e,D_w_e,D_z_c_e,D_z_u_e,D_z_w_e,X_c_e,X_u_e,X_w_e,Kki,Ei,u,v,C,K] = chris_exp_spectral_transport(need_dpxdx)

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

etai = -1:0.2:1;
%Ki = 0.25*[0:10:2*10^3];
Ki = 0.2:0.2:6;

[Kki,Ei] = meshgrid(Ki,etai);

Ei = flipud(Ei);

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
     
     [C,eta,x,delta] = chris_C(C);
    
    delta = delta/dpxdx;
    
    t0 = h_0./hd_0;
    
    tau(i) = (1/1000)*t_rs/t0;
    
    [xm,etam] = meshgrid(x,eta);
    
    area_a =  imgaussfilt(medfilt2(v));
    area_b =  imgaussfilt(medfilt2(v));
    area_a(etam>0)=nan;
    area_b(etam<0)=nan;
    
    area_a(etam<-0.5)=nan;
    area_b(etam>0.5)=nan;
    
    dv = hd_0/(0.28*A_rs);
    
    deltadot = nanmean(max(area_a'))-nanmean(min(area_b'));
    
    if mod(nx,2)==1
        nx=nx-1;
    end
    
    u=imgaussfilt(u(:,1:nx),3)/deltadot;
    v=imgaussfilt(v(:,1:nx),3)/deltadot;
    C=imgaussfilt(C(:,1:nx),3);
    
    h = delta/dpxdx;
    Reh = h*deltadot/nu_m;
    
    de = abs(eta(2)-eta(1));
    
    [K,Z,Pi_w,Pi_u,Pi_c,Pi_z_w,Pi_z_u,Pi_z_c,P_c,P_u,P_w,E_u,E_w,E_c,D_c,D_u,D_z_w,D_z_c,D_z_u,D_w,X_c,X_u,X_w] = chris_exp_spectral_complete([fliplr(u) u fliplr(u)],[fliplr(v) v fliplr(v)],[fliplr(C) C fliplr(C)],h/0.25,dpxdx,de,Reh,0.1);
    
    k = 2*dpxdx*[0:floor(3*nx/2)-1]/(3*nx);
    z = [1:ny]/dpxdx;
    
    [K,E] = meshgrid(log(k),eta);
    K(isinf(K))=0;
    
    filtsize = 2;
    
    E_c = E_c(:,1:floor(3*nx/2));
    E_c_e(:,:,i) = interp2(K,E,E_c,Kki,Ei);
    
    E_u =  E_u(:,1:floor(3*nx/2)) ;
    E_u_e(:,:,i) = interp2(K,E,E_u,Kki,Ei);
    
    E_w =  E_w(:,1:floor(3*nx/2)) ;
    E_w_e(:,:,i) = interp2(K,E,E_w,Kki,Ei);
    
    Pi_c =  Pi_c(:,1:floor(3*nx/2)) ;
    Pi_c_e(:,:,i) = interp2(K,E,Pi_c,Kki,Ei);
    
    Pi_u =  Pi_u(:,1:floor(3*nx/2)) ;
    Pi_u_e(:,:,i) = interp2(K,E,Pi_u,Kki,Ei);
    
    Pi_w =  Pi_w(:,1:floor(3*nx/2)) ;
    Pi_w_e(:,:,i) = interp2(K,E,Pi_w,Kki,Ei);
    
    Pi_z_c =  Pi_z_c(:,1:floor(3*nx/2)) ;
    Pi_z_c_e(:,:,i) = interp2(K,E,Pi_z_c,Kki,Ei);
    
    Pi_z_u =  Pi_z_u(:,1:floor(3*nx/2)) ;
    Pi_z_u_e(:,:,i) = interp2(K,E,Pi_z_u,Kki,Ei);
    
    Pi_z_w =  Pi_z_w(:,1:floor(3*nx/2)) ;
    Pi_z_w_e(:,:,i) = interp2(K,E,Pi_z_w,Kki,Ei);
    
    P_c =  P_c(:,1:floor(3*nx/2)) ;
    P_c_e(:,:,i) = interp2(K,E,P_c,Kki,Ei);
    
    P_u =  P_u(:,1:floor(3*nx/2)) ;
    P_u_e(:,:,i) = interp2(K,E,P_u,Kki,Ei);
    
    P_w =  P_w(:,1:floor(3*nx/2)) ;
    P_w_e(:,:,i) = interp2(K,E,P_w,Kki,Ei);
    
    D_c =  D_c(:,1:floor(3*nx/2)) ;
    D_c_e(:,:,i) = interp2(K,E,D_c,Kki,Ei);
    
    D_u =  D_u(:,1:floor(3*nx/2)) ;
    D_u_e(:,:,i) = interp2(K,E,D_u,Kki,Ei);
    
    D_w =  D_w(:,1:floor(3*nx/2)) ;
    D_w_e(:,:,i) = interp2(K,E,D_w,Kki,Ei);
    
    D_z_c =  D_z_c(:,1:floor(3*nx/2)) ;
    D_z_c_e(:,:,i) = interp2(K,E,D_z_c,Kki,Ei);
    
    D_z_u =  D_z_u(:,1:floor(3*nx/2)) ;
    D_z_u_e(:,:,i) = interp2(K,E,D_z_u,Kki,Ei);
    
    D_z_w =  D_z_w(:,1:floor(3*nx/2)) ;
    D_z_w_e(:,:,i) = interp2(K,E,D_z_w,Kki,Ei);
    
    X_c =  X_c(:,1:floor(3*nx/2)) ;
    X_c_e(:,:,i) = interp2(K,E,X_c,Kki,Ei);
    
    X_u =  X_u(:,1:floor(3*nx/2)) ;
    X_u_e(:,:,i) = interp2(K,E,X_u,Kki,Ei);
    
    X_w =  X_w(:,1:floor(3*nx/2)) ;
    X_w_e(:,:,i) = interp2(K,E,X_w,Kki,Ei);
    
end

end

