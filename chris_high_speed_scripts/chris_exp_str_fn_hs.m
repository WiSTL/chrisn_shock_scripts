function [rx,Sxm,Sx,zeta,alpha] = chris_exp_str_fn_hs(C,n,dpxdx,eta,tau,tau_i)

eta_i = -1:0.01:1;
tau_i = tau_i;

[C_i] = chris_interp_hs(C,eta,tau,eta_i,tau_i);

C_i(isnan(C_i))=0;

[ne,nx,nt] = size(C_i);

ft = fittype('a+b*x');

itop = find(eta_i>-0.5);
range_cp = find(eta_i(itop)<0.5);
range_cp = itop(1):itop(range_cp(end));

for i = 1:nt
     Ci = squeeze(C_i(:,:,i));
     [Ci,eta,x,delta,y0] = chris_C(Ci);
     [rho,mu] = chris_rho(Ci);

    [xm,etam] = meshgrid(x,eta);
    
    Cm = mean(Ci');
    rm = mean(rho');
    
    C_p = Ci;
    r_p = rho;

    
    for j = 1:nx
        C_p(:,j) = Ci(:,j) - Cm';
    end

    [Sax,rax] = chris_no_structure_fn(C_p(range_cp,:),n,dpxdx);
    
    Sx(:,:,i) = Sax;
    Sxm(:,i) = mean(Sax);
    rx(:,i) = rax;
    
    for j = 1:length(range_cp)
        a = log(Sx(j,2:10,i))';  
        a(isnan(a)) = 0;
        a(isinf(a)) = 0;
        f = fit(log(rx(2:10,i)),a, ft,'StartPoint',[1,1]);
        zeta(j,i) = f.b;
        alpha(j,i) = exp(f.a);
    end
    
end

end

