function [Pn,Tn,Dn,en,eta_m] = chris_exp_energy_balance(need_dpxdx,eta_range)

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

%     u = medfilt2(medfilt2(u,[3 3]));
%     v = medfilt2(medfilt2(v,[3 3]));
     [ny,nx] = size(C);
     [C,eta,x,delta,y0] = chris_C(C);
     [rho,mu] = chris_rho(C);

    [xm,etam] = meshgrid(x,eta);

    Cmm = 0.5*erfc(etam);
    Cpm = C-Cmm;
    
    area_a = medfilt2(v);
    area_b = medfilt2(v);
    area_a(etam>0)=nan;
    area_b(etam<0)=nan;
    deltadot = nanmean(max(area_a'))-nanmean(min(area_b'));

    [P,T,D,e,dd,term] = fluc_equation_terms(u,v,C,rho,deltadot,delta/dpxdx,mu./rho,eta);
        
    P12 = interp1(eta,P(:,1),eta_m,'pchip',0);
    T12 = interp1(eta,T(:,1),eta_m,'pchip',0);
    D12 = interp1(eta,D(:,1),eta_m,'pchip',0);
    e12 = interp1(eta,e(:,1),eta_m,'pchip',0);
    
    P22 = interp1(eta,P(:,2),eta_m,'pchip',0);
    T22 = interp1(eta,T(:,2),eta_m,'pchip',0);
    D22 = interp1(eta,D(:,2),eta_m,'pchip',0);
    e22 = interp1(eta,e(:,2),eta_m,'pchip',0);
    
    P32 = interp1(eta,P(:,3),eta_m,'pchip',0);
    T32 = interp1(eta,T(:,3),eta_m,'pchip',0);
    D32 = interp1(eta,D(:,3),eta_m,'pchip',0);
    e32 = interp1(eta,e(:,3),eta_m,'pchip',0);
    
    
    Pn(1,i) = trapz(eta_m,P12);
    Tn(1,i) = trapz(eta_m,T12);
    Dn(1,i) = trapz(eta_m,D12);
    en(1,i) = trapz(eta_m,e12);
    
    Pn(2,i) = trapz(eta_m,P22);
    Tn(2,i) = trapz(eta_m,T22);
    Dn(2,i) = trapz(eta_m,D22);
    en(2,i) = trapz(eta_m,e22);
    
    Pn(3,i) = trapz(eta_m,P32);
    Tn(3,i) = trapz(eta_m,T32);
    Dn(3,i) = trapz(eta_m,D32);
    en(3,i) = trapz(eta_m,e32);
    
end
end

