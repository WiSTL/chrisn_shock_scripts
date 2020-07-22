function [delta,y0] = chris_exp_mix_width_hs(C)

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
    for j = 1:nx
        Ci = squeeze(C(:,:,i));
     
        y = [1:ny];
        x = [1:nx];
    
        ft = fittype('0.5*erfc((4/sqrt(2*pi))*(x-b)/a)');
        f = fit(y(:),Ci(:,j), ft,'StartPoint',[200,90]);

        delta(i,j) = f.a;
        y0(i,j) = f.b;
    end
    
end

end

