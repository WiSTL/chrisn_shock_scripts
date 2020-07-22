function [cm,um,vm,Pi_i,Pi_h,Rm,Etam] = chris_exp_1D_structure_fn_terms(n,need_dpxdx)

[filename,pathname] = uigetfile('*.mat','multiselect','on');

fileno = size(filename);
fileno = fileno(2)

if fileno>1
    for i = 1:fileno
        fileLocations{i} = [pathname filename{i}];
    end 
else
    fileLocations = [pathname filename];
end

r_m = [0:0.0005:0.05];
eta_m = -0.5:0.01:0.5;

[Rm,Etam] = meshgrid(r_m,eta_m);

um = 0*Rm;
vm = 0*Rm;
cm = 0*Rm;
Pi_i = 0*Rm;
Pi_h = 0*Rm;

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
     [delta_dot,v,~] = chris_v(eta,v,C);
     [rho,mu] = chris_rho(C);
     
     u = u/delta_dot;
     v = (v-mean2(v))/delta_dot;
    
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
    

    [Sux,rx,du,ua] = chris_no_structure_fn(u_p,n,dpxdx);
    
    [R,Eta] = meshgrid(rx,eta(1:end));
    
    u_m =  interp2(R,Eta,Sux,Rm,Etam,'cubic',0);
    um = um+imgaussfilt(u_m);
    
    [Svx,rx,dv,va] = chris_no_structure_fn(v_p,n,dpxdx);
    
    [R,Eta] = meshgrid(rx,eta(1:end));
    
    v_m =  interp2(R,Eta,Svx,Rm,Etam,'cubic',0);
    vm = vm+imgaussfilt(v_m);
    
    [Scx,rx,dc,ca] = chris_no_structure_fn(c_p,n,dpxdx);
    
    [R,Eta] = meshgrid(rx,eta(1:end));
    
    c_m =  interp2(R,Eta,Scx,Rm,Etam,'cubic',0);
    cm = cm+imgaussfilt(c_m);
    
    [Swcx,rx,dwc,wca] = chris_no_structure_fn(v_p.*c_p,n,dpxdx);
    
    [Sucx,rx,duc,uca] = chris_no_structure_fn(u_p.*c_p,n,dpxdx);
    
    Pi_inhom = interp2(R,Eta,n*squeeze(mean(dwc.*dc.^(n-1),2)),Rm,Etam,'cubic',0);
    Pi_hom = interp2(R,Eta,n*squeeze(mean(uca.*dc.^(n-1),2)),Rm,Etam,'cubic',0);
    
    Pi_i = Pi_i + imgaussfilt(Pi_inhom);
    Pi_h = Pi_h + imgaussfilt(Pi_hom);
    

end

um = um/fileno;
vm = vm/fileno;
cm = cm/fileno;
Pi_h = Pi_h/fileno;
Pi_i = Pi_i/fileno;

figure(1)
pcolor(Rm,Etam,um)
shading flat
set(gca,'xscale','log');


figure(2)
pcolor(Rm,Etam,vm)
shading flat
set(gca,'xscale','log');


figure(3)
pcolor(Rm,Etam,cm)
shading flat
set(gca,'xscale','log');


figure(4)
pcolor(Rm,Etam,Pi_i)
shading flat
set(gca,'xscale','log');



end

