function [kd,h,cp2,vp2,dv,F,tau] = chris_exp_wave_interaction(need_dpxdx)

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

k_m = [0:0.5:1.2*10^3];
eta_m = -0.5:0.01:0.5;

[Km,Etam] = meshgrid(k_m,eta_m);

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
    vm = mean(v');
    um = mean(u');
    
    y = [1:ny];
    x = [1:nx];

    ft = fittype('0.5*erfc((4/sqrt(2*pi))*(x-b)/a)');
    f = fit(y(:),Cm(:), ft,'StartPoint',[200,90]);

    delta = f.a;
    y0 = f.b;

    eta = (4/sqrt(2*pi))*(y-y0)/delta;
    
    [rho,mu] = chris_rho(C);
    
    rm = mean(rho');
    
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
   
   x = [0:nx-1]/dpxdx;
   
   r_p = C;
   c_p = C;
   v_p = v;
    
    for j = 1:nx
        r_p(:,j) = rho(:,j) - rm';
        c_p(:,j) = C(:,j) - Cm';
        v_p(:,j) = v(:,j) - vm';
    end
    
    for j = 1:ny
    [e(j,:),k] = pspectrum(r_p(j,:),x);
    end
    
    [K,Eta] = meshgrid(k,eta(1:end));
    r_m =  interp2(K,Eta,e,Km,Etam,'cubic',NaN);

    [~,ik] = findpeaks(nanmean(r_m),Km(1,:));
    
    R = (1+A_s)/(1-A_s);
    
    kd(i) = ik(1);
    h(i) = delta/dpxdx;
    cp2(i) = sqrt(mean2(c_p.^2));
    vp2(i) = sqrt(mean2(v_p.^2));
    F(i) = energy_deposited(kd(i)*h_0,R);
    t0 = h_0./hd_0;
    tau(i) = (1/1000)*t_rs/t0;
    dv(i) = hd_0/(0.28*A_rs);
    
    clear e;
end




end

