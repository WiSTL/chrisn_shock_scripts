function [t11,t12,t22,u_f,v_f,K,e,Ret,Ci,csaa,csab,csbb,ce,q2] = LES_profile(u,v,C,rho,mu,dpxdx)

    [ny,nx] = size(C);
    df=nx;

    t11 = mean((rho.*u.^2)') - (mean((rho.*u)').^2)./(mean(rho'));
    t12 = mean((rho.*u.*v)') - (mean((rho.*u)').*mean((rho.*v)'))./(mean(rho'));
    t22 = mean((rho.*v.^2)') - (mean((rho.*v)').^2)./(mean(rho'));
    
    q1 = mean((rho.*u.*C)') - (mean((rho.*C)').*mean((rho.*u)'))./(mean(rho'));
    q2 = mean((rho.*v.*C)') - (mean((rho.*C)').*mean((rho.*v)'))./(mean(rho'));
    
    j1 = mean((rho.*u.*(u.^2 + v.^2))') - (mean((rho.*u)').*mean((rho.*(u.^2 + v.^2))'))./(mean(rho'));
    j2 = mean((rho.*v.*(u.^2 + v.^2))') - (mean((rho.*v)').*mean((rho.*(u.^2 + v.^2))'))./(mean(rho'));
    
    [dudy,dudx] = chris_gradient(imgaussfilt(u,10),1/dpxdx,1/dpxdx);
    [dvdy,dvdx] = chris_gradient(imgaussfilt(v,10),1/dpxdx,1/dpxdx);
    
    S11 = dudx;
    S12 = 0.5*(dudy+dvdx);
    S22 = dvdy;
    
    sigma11 = 2*mu.*S11 - (2/3)*mu.*(S11+S22);
    sigma12 = 2*mu.*S12 - (2/3)*mu.*(S11+S22);
    sigma22 = 2*mu.*S22 - (2/3)*mu.*(S11+S22);
     
    mu_f = mean((rho.*mu)')./mean(rho');
    u_f = mean((rho.*u)')./mean(rho');
    v_f = mean((rho.*v)')./mean(rho');
    rho_f = mean(rho');
    
    S11_f = mean((rho.*S11)')./mean(rho');
    S12_f = mean((rho.*S12)')./mean(rho');
    S22_f = mean((rho.*S22)')./mean(rho');
    
    sigma11_f = 2*mu_f.*S11_f - (2/3)*mu_f.*(S11_f+S22_f);
    sigma12_f = 2*mu_f.*S12_f - (2/3)*mu_f.*(S11_f+S22_f);
    sigma22_f = 2*mu_f.*S22_f - (2/3)*mu_f.*(S11_f+S22_f);
    
    e_v = mean((sigma11.*S11+2*sigma12.*S12+sigma22.*S22)') - sigma11_f.*S11_f+2*sigma12_f.*S12_f+sigma22_f.*S22_f;
    
    C_f = mean((rho.*C)')./mean(rho');
    
    Kf = 0.5*(mean(rho_f.*(u_f.^2 + v_f.^2))/(mean(rho_f)) - (mean(rho_f.*u_f)^2 + mean(rho_f.*v_f)^2)/(mean(rho_f)^2));
    kf = mean(t11+t22)/(2*mean(rho_f));
    
    K = Kf+kf;
    
    t11p = t11-mean(t11);
    t12p = t12-mean(t12);
    t22p = t22-mean(t22);
    
    S11p = S11_f - mean(S11_f);
    S12p = S12_f - mean(S12_f);
    S22p = S22_f - mean(S22_f);
    
    Sigma11p = sigma11_f - mean(sigma11_f);
    Sigma12p = sigma12_f - mean(sigma12_f);
    Sigma22p = sigma22_f - mean(sigma22_f);
    
    eres = mean(Sigma11p.*S11p + Sigma12p.*S12p + Sigma22p.*S22p)/mean(rho_f);
    esgs = -1*mean(t11p.*S11p + t12p.*S12p + t22p.*S22p)/mean(rho_f);
    
    e = eres+esgs;
    
    up = sqrt(2*K/3);
    
    l = (up^3)/e;
    
    Ret = up*l*mean(rho_f)/mean(mu_f);
    
    S_f_abs = sqrt(2*(S11_f.*S11_f+2*(S12_f.*S12_f)+S22_f.*S22_f));
    tau_kk = t11+t22;
    
    Ci = tau_kk./(2*rho_f.*((1/12)*(df/dpxdx)^2).*S_f_abs.^2);

    tau_aa = t11 - tau_kk/3;
    tau_ab = t12 - tau_kk/3;
    tau_bb = t22 - tau_kk/3;

    sfkk = S11_f+S22_f;
    
    csaa = -tau_aa./(2*rho_f.*((1/12)*(df/dpxdx)^2).*S_f_abs.*(S11_f - sfkk/3));
    csab = -tau_ab./(2*rho_f.*((1/12)*(df/dpxdx)^2).*S_f_abs.*(S12_f - sfkk/3));
    csbb = -tau_bb./(2*rho_f.*((1/12)*(df/dpxdx)^2).*S_f_abs.*(S22_f - sfkk/3));

    df2 = ((1/12)*(df/dpxdx)^2);

    ce = e_v./((1/sqrt(df2)).*rho_f.*(df2.*S_f_abs.^2).^3);
end

