function [C,eta,x,delta,y0,delta2,y02] = chris_C(C)
    
    %%%%%%%%%%%%%%%%%%%%%
    %
    %   Input:
    %   C - 2d slice of the mole fracion field
    %   
    %   Output:
    %   C - clipped output at 0 and 1
    %   eta - (z-z0)/delta
    %   x - x coordinate in px
    %   delta - mixing width in px (curve fit)
    %   y0 - mixing layer centroid (curve fit)
    %   delta2 - mixing width in px (integral)
    %   y02 - mixing layer centroid (integral)
    %
    %%%%%%%%%%%%%%%%%%%%%

    C(C>1)=1;
    C(C<0)=0;
     
    [ny,nx] = size(C);
     
    Cm = smooth(squeeze(nanmean((C)')));
    y = [1:ny];
    x = [1:nx];
    
    ft = fittype('0.5*erfc((4/sqrt(2*pi))*(x-b)/a)'); %x-b
    f = fit(y(:),Cm(:), ft,'StartPoint',[200,90]);
    
    delta = f.a;
    y0 = f.b;
    
    
    delta2 = 4*trapz(Cm.*(1-Cm));
    y02 = 4*trapz(Cm.*(1-Cm).*y')/delta;
    
    eta = -1*(y-y0)/delta; %(4/sqrt(2*pi))*

end

