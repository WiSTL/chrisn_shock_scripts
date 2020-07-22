function [Xi] = C_chemical_product(C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Weber taken from Cook and Dimotakis
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = mixing_layer_mask(C,subplus(1-C));

Cm = 2*min(C,1-C);

Cav = nanmean(C');
Cmav = nanmean(Cm');

Cavm = 2*min(Cav,1-Cav);

Xi = nansum(Cmav)/nansum(Cavm);


end

