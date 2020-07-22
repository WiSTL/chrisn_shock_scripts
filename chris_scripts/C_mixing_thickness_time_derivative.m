function [h_dot] = C_mixing_thickness_time_derivative(C,u,v)

meanC = mean(C');

[dcdx,dcdy] = chris_gradient(medfilt2(C),1,1);


meanC_dot = mean((u.*dcdx+v.*dcdy)');


h_dot = 4*sum(meanC_dot.*(2*meanC-1));

end

