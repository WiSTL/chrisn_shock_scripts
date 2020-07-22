function [delta_c,delta_v,delta_av] = chris_fluc_length_scales(eta,cp2,up2,vp2,apvp)

    ft = fittype('a*exp(-b*(x-c)^2)');
    ft1 = fittype('a*(x^2)*exp(-b*(x-c)^2)');
    
    fc = fit(eta(:),cp2(:), ft,'StartPoint',[0.1,1,0]);
    %fu = fit(eta(:),up2(:), ft1,'StartPoint',[0.1,1,0]);
%     fv = fit(eta(:),vp2(:), ft,'StartPoint',[0.1,1,0]);
%     fav = fit(eta(:),apvp(:), ft,'StartPoint',[0.1,1,0]);

    delta_c = (fc.a)*sqrt(pi)/(fc.b);
    %delta_u = (fu.a)*sqrt(pi)*(1+2*(fu.b)*(fu.c)^2)/(fu.b^(3/2));
    delta_v = 0;%(fv.a)*sqrt(pi)/(fv.b);
    delta_av = 0;%(fav.a)*sqrt(pi)/(fav.b);

end

