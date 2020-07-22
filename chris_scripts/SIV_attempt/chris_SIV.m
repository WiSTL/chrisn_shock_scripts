function [u,w] = chris_SIV(C,Z,X,T,V0)

D = 0.00059;

[nz,nx,nt] = size(C);

u = 0*C;
w = 0*C;

dz = abs(squeeze(Z(2,1,1) - Z(1,1,1)))

[dcz,dcx,dct] = chris_gradient_hs(C,X,Z,T);

[d2ctz,d2ctx,~] = chris_gradient_hs(dct,X,Z,T);

[d2cxz,d2cx2,d2cxt] = chris_gradient_hs(dcx,X,Z,T);

[d2cz2,d2czx,d2czt] = chris_gradient_hs(dcz,X,Z,T);

d2czt = 0.5*medfilt3(d2czt+d2ctz);
d2cxt = 0.5*medfilt3(d2cxt+d2ctx);
d2cxz = 0.5*medfilt3(d2cxz+d2czx);
d2cz2 = medfilt3(d2cz2);
d2cx2 = medfilt3(d2cx2);

for n = 1:nt
    u(:,:,n) = V0(n);
end


for n = 1:nt
    for i = 1:nx
      for j = 2:nz
         u(j,i,n) = dz*((d2cxz(j-1,i,n))*u(j-1,i,n) + (d2cz2(j-1,i,n))*w(j-1,i,n) + d2czt(j-1,i,n)); 
         
         w(j,i,n) = dz*((-d2cx2(j-1,i,n))*u(j-1,i,n) + (d2cxz(j-1,i,n))*w(j-1,i,n) - d2cxt(j-1,i,n));
      end
    end
end

end

