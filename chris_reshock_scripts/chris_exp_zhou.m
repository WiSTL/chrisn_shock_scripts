function [mass,time,time1,tau,tau1] = chris_exp_zhou(need_dpxdx)

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

mass=[];
time=[];
tau=[];
 
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
    
    M0 = 14.8;
    M1 = 40;
    
    W = C*M0./(C*M0+(1-C)*M1);
    
    Cmm = 0.5*erfc(-etam);
    Cpm = C-Cmm;
    
    deltadot = 1;%delta_dot;
    
    dx = abs(x(1) - x(2))/dpxdx;
    deta = abs(eta(1)-eta(2));
    delta = delta/dpxdx;
    
    a = deta*trapz(W.*(1-W).*rho);
    b = dx*trapz(a);
    
    mass(i) = delta*b;
    
    time(i) = (delta/delta_dot)*((100/5)^2)*(1.6*10^5)^-0.5;
    time1(i) = (h_0/hd_0)*((100/5)^2)*(1.6*10^5)^-0.5;
    
    tau(i) = (1/1000)*t_rs/time(i);

    
    tau1(i) = (1/1000)*t_rs/time1(i);
    
end
end

