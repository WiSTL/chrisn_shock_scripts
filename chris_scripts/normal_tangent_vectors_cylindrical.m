function [Tr,Tt,Nr,Nt] = normal_tangent_vectors_cylindrical(theta,r)

%%%%%%%%%%%%%%%%%%%
%%% Tangent
%%%%%%%%%%%%%%%%%%%

Tr = diff(r)./diff(theta);
Tt = r(2:end);

Tmag = sqrt(Tr.^2 + Tt.^2);

Tr = Tr./Tmag;
Tt = Tt./Tmag;

%%%%%%%%%%%%%%%%%%%
%%% normal
%%%%%%%%%%%%%%%%%%%

Nr = diff(Tr)./diff(theta(2:end)) - Tt(2:end);
Nt = Tr(2:end)+diff(Tt)./diff(theta(2:end));

Nmag = sqrt(Nr.^2 + Nt.^2);

Nr = Nr./Nmag;
Nt = Nt./Nmag;


end

