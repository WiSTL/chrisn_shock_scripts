function [Mt,R,cp20,Mi] = chris_exp_theory_energy_deposition(need_dpxdx)

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

Mt = [];
R = [];

for i = 1:fileno
    i
    
     if fileno>1
        file = fileLocations{i};
    else
        file = fileLocations;
     end
   
     load(file);
     
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
    
    vm = mean2(v);
    v = v-vm;
    u = u;
    
    Cm = mean(C');
    rm = mean(rho');
    
    C_p = C;
    v_p = v;
    u_p = u;
    r_p = rho;
    vm2 = mean(v');
    um2 = mean(u');
    
    for j = 1:nx
        C_p(:,j) = C(:,j) - Cm';
    end
    
    for j = 1:nx
        r_p(:,j) = rho(:,j) - rm';
    end
    
    for j = 1:nx
        v_p(:,j) = v(:,j) - vm2';
    end
    
    for j = 1:nx
        u_p(:,j) = v(:,j) - um2';
    end
    
    Mt(i) = max(sqrt(max((v_p.^2)')))/450;
    R(i) = (1+A_rs)/(1-A_rs);
    cp20(i) = max(sqrt(mean((C_p.^2)')));
    Mi(i) = M_rs;
    
end


end

