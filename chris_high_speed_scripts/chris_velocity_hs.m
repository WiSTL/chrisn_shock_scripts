function [u,w,err] = chris_velocity_hs(C,X,Z,T,V0)

[dcz,dcx,dct] = chris_gradient_hs(C,X,Z,T);

dcz = medfilt3(dcz);
dcx = medfilt3(dcx);
dct = medfilt3(dct);

[nz,~,nt] = size(C);

dz = Z(2,1,1)-Z(1,1,1);

ug = 0*C;
for i=1:nt
    wg(:,:,i) = 0*C(:,:,i) + V0(i);
end

w1 = wg;
w2 = wg;

for k = 1:20
    k
    u = -1*(wg.*dcz+dct)./dcx;
    u(isnan(u))=0;
    u = medfilt3(u);
    
    w = -1*(u.*dcx+dct)./dcz;
    w(isnan(w))=0;
    w = medfilt3(w);
    
    [~,dux,~] = chris_gradient_hs(u,X,Z,T);
    for i = 2:nz
        w1(i,:,:) = wg(i-1,:,:) - dz*dux(i,:,:);
        w2(end+2-i,:,:) = wg(end+1-i,:,:) - dz*dux(end+2-i,:,:);
    end

    w1(isnan(w1))=0;
    w2(isnan(w2))=0;

    w = 0.5*(w1+w2);

    err(k) = nanmean(dct(:) + u(:).*dcx(:) + w(:).*dcz(:));
    
    ug = u;
    wg = w;

end

end

