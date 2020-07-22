function [x1,A,lcx,lcy,eps_c] = chris_exp_pdfs_hs(C)

[ny,nx,nt] = size(C);

eta_m = -0.5:0.01:0.5;

xm1 = 0.05:0.01:0.95;

A = zeros(length(xm1),nt);
x1 = zeros(length(xm1),nt);

for i = 1:nt
     Ci = squeeze(C(:,:,i));
     [Ci,eta,x,delta,y0] = chris_C(Ci);
     [rho,mu] = chris_rho(Ci);

    [xm,etam] = meshgrid(x,eta);
    
    [h,i5,i95] = concentration_centerline_5_95(Ci);
    
    Cm = mean(Ci');
    rm = mean(rho');
    
    C_p = Ci;
    r_p = rho;

    
    for j = 1:nx
        C_p(:,j) = Ci(:,j) - Cm';
    end
    
    for j = 1:nx
        r_p(:,j) = rho(:,j) - rm';
    end
    
    [dcpdy,dcpdx] = chris_gradient(C_p,1,1);
    
    lcx(:,i) = sqrt(mean((C_p.^2)')./(mean((dcpdx.^2)')));
    lcy(:,i) = sqrt(mean((C_p.^2)')./(mean((dcpdy.^2)')));  
    
    eps_c(:,i) = mean((dcpdx.^2+dcpdy.^2)');
    
    itop = find(eta>-1);
    range_cp = find(eta(itop)<1);
    
    [A1,x2] = chris_pdf(Ci(range_cp,:),0.01,0.05,0.95);
    
    Am = interp1(x2,A1,xm1,'pchip',0);
    
    A(:,i) = Am;
    x1(:,i) = xm1;
    
    
    
end

end

