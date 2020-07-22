function [f_x,dc_x,f_y,dc_y,A_x,A_y,Csisx,Csisy] = scalar_increment_statistics(C,r,dc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%      Scalar Increment Statistics - x direction
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = round(C/max(max(C)),3);

[nx,ny] = size(C);

A_x=[];
n=0;

I2 = zeros(nx+r,ny);
I3 = zeros(nx+r,ny);

for i=1+r:nx+r
    for j=1:ny
        I2(i,j)=C(i-r,j);
        I3(i-r,j)=C(i-r,j);
    end
end

Csisx = I2-I3;

dcx = dc;

for i=dcx
    n=n+1;
    [A_x(n),~] = size(Csisx(Csisx==i));    
end


for i=1:n
   if(A_x(i)==0)
       A_x(i) = nan;
       dcx(i) = nan;
   end
end

A_x=A_x(~any(isnan(A_x),1));%Area for a given Concentration
dc_x=dcx(~any(isnan(dcx),1));

A_x = interp1(dc_x,A_x,[-0.15:0.0005:0.15]);

dc = [-0.15:0.0005:0.15];
dc_x = dc;

ddc = diff(dc);
dA = diff(smooth(A_x))';

f_x = abs((dA)./ddc/(nx*ny));%Probability Density Function

% f_x=log(f_x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%      Scalar Increment Statistics - y direction
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_y=[];
n=0;

I2 = zeros(nx,ny+r);
I3 = zeros(nx,ny+r);

for i=1:nx
    for j=1+r:ny+r
        I2(i,j)=C(i,j-r);
        I3(i,j-r)=C(i,j-r);
    end
end

Csisy = I2-I3;

dcy = dc;

for i=dcy
    n=n+1;
    [A_y(n),~] = size(Csisy(Csisy==i));    
end

for i=1:n
   if(A_y(i)==0)
       A_y(i) = nan;
       dcy(i) = nan;
   end
end

A_y=A_y(~any(isnan(A_y)));%Area for a given Concentration
dc_y=dcy(~any(isnan(dcy)));

ddc = diff(dc_y);
dA = diff(smooth(A_y),25);

f_y = 0;

% f_y = 0;abs((dA)./ddc/(nx*ny));%Probability Density Function

% f_y=log(f_y);


end

