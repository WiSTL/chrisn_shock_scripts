function [d,ds] = chris_exp_length_scales_hs(C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   1 - mixing thickness
%   2 - displacment thickness
%   3 - Momentum Thickness
%   4 - energy thickness
%   5 - stress thickness, ds - h_dot
%   6 - dissipation thickness
%   7-9 - taylor scale
%
%   d = integral <> - <>^2 dz
%   ds = integral < - ^2> dz - d
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ny,nx,nt] = size(C);

eta_m = -0.5:0.01:0.5;

for i = 1:nt     
    Ci = squeeze(C(:,:,i));
     [~,eta,x,delta,y0] = chris_C(Ci);
     [rho,mu] = chris_rho(Ci);
     
    [xm,etam] = meshgrid(x,eta);

    Cmm = 0.5*erfc(etam);
    Cpm = Ci-Cmm;
    rm = mean(rho');
    
    Cm = mean(Ci');
    
    d12 = interp1(eta,mean((Ci.*(1-Ci))'),eta_m,'pchip',0);

    d(1,i) = delta;
    
    d(2,i) = delta-delta*trapz(eta_m,d12);

    C_p = Ci;
    
    for j = 1:nx
        C_p(:,j) = Ci(:,j) - Cm';
    end
    
    itop = find(eta>-0.5);
    range_cp = find(eta(itop)<0.5);
    
    if length(range_cp)>4
    [~,~,~,~,l_cx,l_cy] = chris_taylor_microscale(C_p(range_cp,:),C_p(range_cp,:),C_p(range_cp,:));
%     lcxn(i) = l_cx;
%     lcyn(i) = l_cy;

    ds(1,i) = l_cx;
    ds(2,i) = l_cy;
    end

    
end

end

