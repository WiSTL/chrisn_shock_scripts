function [ui,wi,uo,wo,un,wn,vor,dcz,dcx,dct,error1] = chris_velocity_hs_8(C)

dpxdx = 1;%2274.375*2;
dt = 1;%5E-5;

[~,~,nt] = size(C);

[dcz,dcx,dct] = chris_gradient_hs(C,1,1,dt);

textprogressbar('calculating outputs: ');
 
for i = 1:nt
    textprogressbar(100*(i-1)/((nt)));
    [u1,w1,~,~,~,~,~,~,error1]=OpticalFlowPhysics_fun_4(squeeze(C(:,:,i)),squeeze(dcx(:,:,i)),squeeze(dcz(:,:,i)),squeeze(dct(:,:,i)),dpxdx,1,1);

%     u1 = u1*96/0.002651;
%     w1 = w1*96/0.002651;
    
    vor(:,:,i) = chris_vort(u1,w1,1/dpxdx);
    
    [u2,w2] = chris_biot_savart(squeeze(vor(:,:,i)),1/dpxdx,1/dpxdx);
    [~,u3,w3,~] = chris_solvePoissonSOR(squeeze(vor(:,:,i)),5000,1*10^-8,squeeze(vor(:,:,i)),1/dpxdx,1/dpxdx);
    
    uo(:,:,i) = u1;
    wo(:,:,i) = w1;
    
    ui(:,:,i) = u2;
    wi(:,:,i) = w2;
    
    un(:,:,i) = u3;
    wn(:,:,i) = w3;
end

textprogressbar('      done');

% for i = 1:nt-2
%     Df(:,:,1) = squeeze(uo(:,:,i+1));
%     Df(:,:,2) = squeeze(wo(:,:,i+1));
%     
%     Cf = imwarp(squeeze(C(:,:,i+1)),Df/2);
%     C1(:,:,i) = Cf(:,:,i);
% %     
% %     Tf = imwarp(squeeze(T(:,:,i)),Df/2);
% %     Tb = imwarp(T(:,:,i+1),Dfb/2);
% %     T1(:,:,i) = (Tf+Tb)/2;
% end

end

