function [Ty,Tx,Ny,Nx] = normal_tangent_vectors_cartesian(x,y)

%%%%%%%%%%%%%%%%%%%
%%% Tangent
%%%%%%%%%%%%%%%%%%%

Ty = diff(y)./diff(x);
Tx = 1;

Tmag = sqrt(Ty.^2 + Tx.^2);

Ty = Ty./Tmag;
Tx = Tx./Tmag;

%%%%%%%%%%%%%%%%%%%
%%% normal
%%%%%%%%%%%%%%%%%%%

Ny = diff(Ty)./diff(x(2:end));
Nx = diff(Tx)./diff(x(2:end));

Nmag = sqrt(Ny.^2 + Nx.^2);

Ny = Ny./Nmag;
Nx = Nx./Nmag;


end

