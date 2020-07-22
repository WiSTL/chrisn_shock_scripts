function [u1,w1,u,w,vor,div,p,C2,u2,w2,vor2,div2] = chris_velocity_hs_6(C)

dpxdx = 2274.375*2;
dt = 5E-5;

[~,~,nt] = size(C);

for i = 1:nt-1
    i
    [u(:,:,i),w(:,:,i),vor(:,:,i),div(:,:,i),u2(:,:,i),w2(:,:,i),vor2(:,:,i),div2(:,:,i)]=OpticalFlowPhysics_fun(squeeze(C(:,:,i)),squeeze(C(:,:,i+1)),200,2000);
    if i>1
    [p(:,:,i),u1(:,:,i),w1(:,:,i)] = chris_solvePoissonSOR(squeeze(vor(:,:,i)),10000,1*10^-8,squeeze(p(:,:,i-1)),1,1);
    else
    [p(:,:,i),u1(:,:,i),w1(:,:,i)] = chris_solvePoissonSOR(squeeze(vor(:,:,i)),10000,1*10^-8,0,1,1);
    end
end

for i = 1:nt-1
    Df(:,:,1) = squeeze(u(:,:,i));
    Df(:,:,2) = squeeze(w(:,:,i));
    
    Cf = imwarp(squeeze(C(:,:,i)),Df/2);
    C2(:,:,i) = Cf;
end

end

