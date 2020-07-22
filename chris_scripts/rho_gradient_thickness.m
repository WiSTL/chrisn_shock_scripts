function [h] = rho_gradient_thickness(r)

[drdx,~] = chris_gradient(medfilt2(r),0.00019,0.00019);

a = abs(drdx)./r;

h = 1/mean(a(:));

end

