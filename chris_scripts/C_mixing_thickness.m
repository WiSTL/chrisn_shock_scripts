function [h] = C_mixing_thickness(C)

meanC = mean(C');

h = 4*sum(meanC.*(1-meanC));

end

