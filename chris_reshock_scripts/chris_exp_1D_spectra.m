function [] = chris_exp_1D_spectra()

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
    
    if i==15 || i==16 || i==17
        dpxdx = 2274.375*2;
    else
        dpxdx = 2274.375;
    end

    r1 = 4.6081;
    r2 = 7.622;
    
    rho = r2+(r1-r2)*C;
    
    figure(1);
    E = u.^2 + v.^2;
    [e,k] = power_spectra_1D(E,1/dpxdx,1,1);
    a = mean(e);
    loglog(k,a/max(a))
    title('Energy Spectra in y direction')
    grid on
    hold on
    xlabel('k')
    ylabel('E_k(k)')
    
    figure(2)
    [e,k] = power_spectra_1D(E',1/dpxdx,1,1);
    a = mean(e);
    loglog(k,a/max(a))
    title('Energy Spectra in x direction')
    grid on
    hold on
    xlabel('k')
    ylabel('E_k(k)')
    
    figure(3);
    [e,k] = power_spectra_1D(C.^2,1/dpxdx,1,1);
    a = mean(e);
    loglog(k,a/max(a))
    title('Scalar power Spectra in y direction')
    grid on
    hold on
    xlabel('k')
    ylabel('E_{\xi}(k)')

    figure(4)
    [e,k] = power_spectra_1D((C.^2)',1/dpxdx,1,1);
    a = mean(e);
    loglog(k,a/max(a))
    grid on
    hold on
    title('Scalar power Spectra in x direction')
    xlabel('k')
    ylabel('E_{\xi}(k)')
end


end

