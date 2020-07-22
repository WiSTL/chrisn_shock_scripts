function [] = chris_exp_str_f()

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
    
    figure(1);
    [S2,r] = chris_2D_structure_fn(u,2);
    plot(r/dpxdx,S2/S2(1))
    title('u structure function')
    grid on
    hold on
    xlabel('r')
    ylabel('S_u_2(r/\delta)')
  
    figure(2);
    [S2,r] = chris_2D_structure_fn(v,2);
    plot(r/dpxdx,S2/S2(1))
    title('v structure function')
    grid on
    hold on
    xlabel('r')
    ylabel('S_v_2(r/\delta)')
    
    figure(3);
    [S2,r] = chris_2D_structure_fn(C,2);
    plot(r/dpxdx,S2/S2(1))
    title('\xi structure function')
    grid on
    hold on
    xlabel('r')
    ylabel('S_{\xi2}(r/\delta)')
end


end

