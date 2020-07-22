function [ui,wi,uo,wo,un,wn,vor] = chris_velocity_hs_7(C)

dpxdx = 2274.375*2;
dt = 5E-5;

[~,~,nt] = size(C);

 textprogressbar('calculating outputs: ');

for i = 2:nt-1
    textprogressbar(100*(i-2)/((nt-1)-2));
    [u1,w1]=OpticalFlowPhysics_fun_3(squeeze(C(:,:,i-1)),squeeze(C(:,:,i)),squeeze(C(:,:,i+1)),10,20);

%     u1 = u1*96/0.002651;
%     w1 = w1*96/0.002651;
    
    vor(:,:,i) = chris_vort(u1,w1);
    
    [u2,w2] = chris_biot_savart(squeeze(vor(:,:,i)),1,1);
    [~,u3,w3,~] = chris_solvePoissonSOR(squeeze(vor(:,:,i)),5000,1*10^-8,squeeze(vor(:,:,i)),1,1);
    
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

