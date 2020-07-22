function [cm,um,vm,Re,eps,scalar_eps,Rm,Etam] = chris_exp_1D_structure_fn(n,need_dpxdx)

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

r_m = [0:0.001:0.2];
eta_m = -0.5:0.01:0.5;

[Rm,Etam] = meshgrid(r_m,eta_m);

um = 0*Rm;
vm = 0*Rm;
cm = 0*Rm;

eps = [];
scalar_eps = [];
Re = [];

for i = 1:fileno
    i
    
     if fileno>1
        file = fileLocations{i};
    else
        file = fileLocations;
     end
   
     load(file);
     
     [ny,nx] = size(C);
     
     [C,eta,x,delta,y0] = chris_C(C);
     [delta_dot,v,vm] = chris_v(eta,v,C);
     [rho,mu] = chris_rho(C);
     
     vm = mean2(v);
     v = v-vm;
     
     Cm = mean(C');
    y = [1:ny];
    x = [1:nx];
    
    if need_dpxdx
        if i==15 || i==16 || i==17
            dpxdx = 2274.375*2;
            dx_les = 10;
        else
            dpxdx = 2274.375;
            dx_les = 5;
        end
    end
    
    cmean = mean(C');
    vmean = mean(v');
    umean = mean(u');
    
    c_p = C;
    v_p = v;
    u_p = u;
    
    for j = 1:nx
        c_p(:,j) = C(:,j) - cmean';
    end
    
    for j = 1:nx
        v_p(:,j) = v(:,j) - vmean';
    end
    
    for j = 1:nx
        u_p(:,j) = u(:,j) - umean';
    end
    
    cp2 = mean(c_p(:).^2);
    
    Re(i) = (delta/dpxdx)*delta_dot./nu_m;
    
    delta_dot = mean2(sqrt(u_p.^2+v_p.^2));
   
    x = [0:nx-1]/dpxdx;
   
    [Sux,rx] = chris_no_structure_fn(u_p,n,dpxdx);
    
    [R,Eta] = meshgrid(rx,eta(1:end));
    
    u_m =  interp2(R,Eta,Sux,Rm,Etam,'cubic',NaN);
    um = um+u_m;
    
    [Svx,rx] = chris_no_structure_fn(v_p,n,dpxdx);
    
    [R,Eta] = meshgrid(rx,eta(1:end));
    
    v_m =  interp2(R,Eta,Svx,Rm,Etam,'cubic',NaN);
    vm = vm+v_m;
    
    eps(i) = nanmean((5/4)*nanmean(v_m(:,2:4),1)./r_m(2:4))*(delta/dpxdx)/delta_dot^3;
    
    
    [Scx,rx] = chris_no_structure_fn(c_p,n,dpxdx);
    
    [R,Eta] = meshgrid(rx,eta(1:end));
    
    c_m =  interp2(R,Eta,Scx,Rm,Etam,'cubic',NaN);
    cm = cm+c_m;
    
    scalar_eps(i) = nanmean(nanmean((5/4)*nanmean(c_m(:,2:4),1)./r_m(2:4)))*(delta/dpxdx)/(delta_dot*cp2);
    
end

um = um/fileno;
vm = vm/fileno;
cm = cm/fileno;

figure(1)
pcolor(Rm,Etam,um)
shading flat
set(gca,'xscale','log');
set(gca,'zscale','log');

figure(2)
pcolor(Rm,Etam,vm)
shading flat
set(gca,'xscale','log');
set(gca,'zscale','log');

figure(3)
pcolor(Rm,Etam,cm)
shading flat
set(gca,'xscale','log');
set(gca,'zscale','log');


end

