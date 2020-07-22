function [cm,em,un,vn,Km,Etam] = chris_exp_1D_spectra_spatial_2018a(need_dpxdx)

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

em = 0*Km;
un = 0*Km;
vn = 0*Km;
cm = 0*Km;
rhom = 0*Km;
ration = 0*Km;

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
   
   c_p = C;
   r_p = rho;
   v_p = v;
   u_p = u;
   
   rm = mean(rho');
   
   cmm=[];
   c_p=[];
   r_p=[];
   v_p=[];
   u_p=[];
   
    for j = 1:nx
        cmm(:,j) = Cm';
        c_p(:,j) = C(:,j) - Cm';
        r_p(:,j) = rho(:,j) - rm';
        v_p(:,j) = v(:,j) - vm';
        u_p(:,j) = u(:,j) - um';
    end
   
   
    for j = 1:ny
    [pu(j,:),k] = pspectrum(sqrt(rho(j,:)).*u_p(j,:),x);
    [pv(j,:),~] = pspectrum(sqrt(rho(j,:)).*v_p(j,:),x);
    [Ev(j,:),~] = pspectrum(v_p(j,:),x);
    [Eu(j,:),~] = pspectrum(u_p(j,:),x);
    [cpvp(j,:),~] = pspectrum(c_p(j,:).*v_p(j,:),x);
    [cpup(j,:),~] = pspectrum(c_p(j,:).*u_p(j,:),x);
    [pr(j,:),~] = pspectrum(r_p(j,:),x);
    e(j,:) = pu(j,:)+pv(j,:);
    end
    
    [K,Eta] = meshgrid(k,eta(1:end));
    
    e_m =  interp2(K,Eta,e,Km,Etam,'cubic',NaN);
    em = em+e_m;
    
    r_m =  interp2(K,Eta,pr,Km,Etam,'cubic',NaN);
    rhom = rhom+r_m;
    
    u_m =  interp2(K,Eta,pu,Km,Etam,'cubic',NaN);
    un = un+u_m;
    
    v_m =  interp2(K,Eta,pv,Km,Etam,'cubic',NaN);
    vn = vn+v_m;
    
    ration = ration + medfilt2(v_m.*r_m);
    
    for j = 1:ny
    [Ec(j,:),k] = pspectrum(c_p(j,:),x);
    end
    
    [K,Eta] = meshgrid(k,eta(1:end));
    c_m =  interp2(K,Eta,Ec,Km,Etam,'cubic',NaN);
    cm = cm+c_m;

    clear Ec;
    clear pu;
    clear pv;
    clear pr;
end

em = em/fileno;
cm = cm/fileno;
un = un/fileno;
rhom = rhom/fileno;
ration = ration/fileno;

figure(1)
surf(Km,Etam,log(real(em)))
shading interp
set(gca,'xscale','log');

figure(2)
surf(Km,Etam,log(real(rhom)))
shading interp
set(gca,'xscale','log');

figure(3)
surf(Km,Etam,real(ration))
shading interp
set(gca,'xscale','log');


end

