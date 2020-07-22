function [] = chris_exp_2D_spectra()

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
    
    dpxdx = 2274.375;
    
    figure(1);
    E = u.^2 + v.^2;
    [~, ~, ~, ~, S2dmirror, kxfull, kyfull] = power_spectra_2D(E,1/dpxdx,1/dpxdx);
    contour(kyfull,kxfull,log(S2dmirror),3),shading flat
    title('2D Energy Spectra')
    grid on
    hold on
    xlabel('k_x(m^{-1})')
    ylabel('k_y(m^{-1})')
    
    figure(3);
    [~, ~, ~, ~, S2dmirror, kxfull, kyfull] = power_spectra_2D(C,1/dpxdx,1/dpxdx);
    contour(kyfull,kxfull,log(S2dmirror),3),shading flat
    title('2D scalar power Spectra')
    grid on
    hold on
    xlabel('k_x(m^{-1})')
    ylabel('k_y(m^{-1})')

end


end

