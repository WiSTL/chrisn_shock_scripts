function [cm,em,Km,Etam] = chris_exp_1D_spectra_spatial(need_dpxdx)

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

k_m = [100:0.5:1.4*10^4];
eta_m = -0.5:0.01:0.5;

[Km,Etam] = meshgrid(k_m,eta_m);

em = 0*Km;
cm = 0*Km;

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
    y = [1:ny];
    x = [1:nx];

    ft = fittype('0.5*erfc((4/sqrt(2*pi))*(x-b)/a)');
    f = fit(y(:),Cm(:), ft,'StartPoint',[200,90]);

    delta = f.a;
    y0 = f.b;

    eta = (4/sqrt(2*pi))*(y-y0)/delta;

    [xm,etam] = meshgrid(x,eta);

    Cmm = 0.5*erfc(etam);
    Cpm = C-Cmm;
    
    if need_dpxdx
        if i==15 || i==16 || i==17
            dpxdx = 2274.375*2;
            dx_les = 10;
        else
            dpxdx = 2274.375;
            dx_les = 5;
        end
    end

    r1 = 4.6081;
    r2 = 7.622;
    
    rho = r2+(r1-r2)*C;
    
    [~,dru2dx] = chris_gradient(rho.*u.*u,1/dpxdx,1/dpxdx);
    
    E = u.^2 + v.^2;
    
%     E = dru2dx;
    
    [e,k] = power_spectra_1D(E,1/dpxdx,1,1);
    [K,Eta] = meshgrid(k,eta(1:end-1));
    e_m =  interp2(K,Eta,e,Km,Etam,'cubic',NaN);
    em = em+e_m;
    
    [e,k] = power_spectra_1D(C.^2,1/dpxdx,1,1);
    [K,Eta] = meshgrid(k,eta(1:end-1));
    c_m =  interp2(K,Eta,e,Km,Etam,'cubic',NaN);
    cm = cm+c_m;

end

em = em/fileno;
cm = cm/fileno;

figure(1)
pcolor(Km,Etam,real(em))
shading flat
set(gca,'xscale','log');

figure(2)
pcolor(Km,Etam,real(cm))
shading flat
set(gca,'xscale','log');


end

