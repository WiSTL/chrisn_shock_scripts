% 
% eta1 = [-2:0.001:2];
% eta2 = [-1.6:0.001:1];
% eta3 = [-1.2:0.001:0.8];
% eta4 = [-0.8:0.001:0.5];
% 
% [L1,M1,N1,Cm1,Cp21,Um1,Up21,Wm1,Wp21,Lambda21] = chris_exp_norm_moments(1,eta1);
% [L2,M2,N2,Cm2,Cp22,Um2,Up22,Wm2,Wp22,Lambda22] = chris_exp_norm_moments(0,eta2);
% [L3,M3,N3,Cm3,Cp23,Um3,Up23,Wm3,Wp23,Lambda23] = chris_exp_norm_moments(0,eta3);
% [L4,M4,N4,Cm4,Cp24,Um4,Up24,Wm4,Wp24,Lambda24] = chris_exp_norm_moments(0,eta4);

figure(4)
chris_plot_mean_std(eta1,(Cm1)','k',200);
hold on
chris_plot_mean_std(eta2,(Cm2)','b',200);
hold on
chris_plot_mean_std(eta3,(Cm3)','g',200);
hold on
chris_plot_mean_std(eta4,(Cm4)','r',200);

figure(5)
chris_plot_mean_std(eta1,(Cp21)','k',200);
hold on
chris_plot_mean_std(eta2,(Cp22)','b',200);
hold on
chris_plot_mean_std(eta3,(Cp23)','g',200);
hold on
chris_plot_mean_std(eta4,(Cp24)','r',200);


figure(6)
chris_plot_mean_std(eta1,(Um1)','k',200);
hold on
chris_plot_mean_std(eta2,(Um2)','b',200);
hold on
chris_plot_mean_std(eta3,(Um3)','g',200);
hold on
chris_plot_mean_std(eta4,(Um4)','r',200);

figure(7)
chris_plot_mean_std(eta1,(Up21)','k',200);
hold on
chris_plot_mean_std(eta2,(Up22)','b',200);
hold on
chris_plot_mean_std(eta3,(Up23)','g',200);
hold on
chris_plot_mean_std(eta4,(Up24)','r',200);

figure(8)
chris_plot_mean_std(eta1,(Wm1)','k',200);
hold on
chris_plot_mean_std(eta2,(Wm2)','b',200);
hold on
chris_plot_mean_std(eta3,(Wm3)','g',200);
hold on
chris_plot_mean_std(eta4,(Wm4)','r',200);

figure(9)
chris_plot_mean_std(eta1,(Wp21)','k',200);
hold on
chris_plot_mean_std(eta2,(Wp22)','b',200);
hold on
chris_plot_mean_std(eta3,(Wp23)','g',200);
hold on
chris_plot_mean_std(eta4,(Wp24)','r',200);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %   L_ij
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % A=[];
% % 
% % for i = 1:5
% %     for j = 1:5
% %       La1 = squeeze(L1(i,j,:,:));
% %       [nex1] = size(La1);
% %       La2 = squeeze(L2(i,j,:,:));
% %       [nex2] = size(La2);
% %       La3 = squeeze(L3(i,j,:,:));
% %       [nex3] = size(La3);
% %       La4 = squeeze(L4(i,j,:,:));
% %       [nex4] = size(La4);
% %       
% %       La = [La1(:)' La2(:)' La3(:)' La4(:)'];
% %       
% %       Ca=[];
% %       for n = 1:nex1
% %           Ca = [Ca Cm1(1,:)];
% %       end
% %       for n = 1:nex2
% %           Ca = [Ca Cm2(1,:)];
% %       end
% %       for n = 1:nex3
% %           Ca = [Ca Cm3(1,:)];
% %       end
% %       for n = 1:nex4
% %           Ca = [Ca Cm4(1,:)];
% %       end
% %       
% %       [A(i,j,:,:),x,y] = chris_jpdf(Ca(:),La(:),20,[0.05 -10],[0.95 10]);
% %       
% %     end
% % end
% 
% 
% figure(1)
% 
% for i = 1:5
%     for j = 1:5
%         subplot(5,5,i+5*(j-1))
%         
%         cma = [nanmean(Cm2,1) nanmean(Cm3,1) nanmean(Cm4,1)];
%         lm = [nanmean(squeeze(L2(i,j,:,:)),1) nanmean(squeeze(L3(i,j,:,:)),1) nanmean(squeeze(L4(i,j,:,:)),1)];
%         
%         plot(cma,lm,'.','MarkerSize',4)
%         
% %         plot(nanmean(Cm1,1),nanmean(squeeze(L1(i,j,:,:)),1),'+')
% %         hold on
% %         plot(nanmean(Cm2,1),nanmean(squeeze(L2(i,j,:,:)),1),'+')
% %         plot(nanmean(Cm3,1),nanmean(squeeze(L3(i,j,:,:)),1),'+')
% %         plot(nanmean(Cm4,1),nanmean(squeeze(L4(i,j,:,:)),1),'+')
%         
%         
%         axis([0 1 -5 5])
%         axis square
%         set(gca,'fontSize',16)
% 
%     end
% end
% 
% subplot(5,5,1)
% ylabel('0','interpreter','latex')
% subplot(5,5,6)
% ylabel('1','interpreter','latex')
% subplot(5,5,11)
% ylabel('2','interpreter','latex')
% subplot(5,5,16)
% ylabel('3','interpreter','latex')
% subplot(5,5,21)
% ylabel('4','interpreter','latex')
% xlabel('0','interpreter','latex')
% subplot(5,5,22)
% xlabel('1','interpreter','latex')
% subplot(5,5,23)
% xlabel('2','interpreter','latex')
% subplot(5,5,24)
% xlabel('3','interpreter','latex')
% subplot(5,5,25)
% xlabel('4','interpreter','latex')
% 
% sgtitle('$L_{i,j}$','interpreter','latex')
% 
% 
% 
figure(2)

for i = 1:7
    for j = 1:7
        subplot(7,7,i+7*(j-1))
        
        cma = [nanmean(Cm2,1) nanmean(Cm3,1) nanmean(Cm4,1)];
        lm = [nanmean(squeeze(M2(i,j,:,:)),1) nanmean(squeeze(M3(i,j,:,:)),1) nanmean(squeeze(M4(i,j,:,:)),1)];
        
        plot(cma,lm,'.','MarkerSize',4)
        
        plot(nanmean(Cm1,1),nanmean(squeeze(M1(i,j,:,:)),1),'k.')
        hold on
        plot(nanmean(Cm2,1),nanmean(squeeze(M2(i,j,:,:)),1),'g.')
        plot(nanmean(Cm3,1),nanmean(squeeze(M3(i,j,:,:)),1),'b.')
        plot(nanmean(Cm4,1),nanmean(squeeze(M4(i,j,:,:)),1),'r.')
        
        axis([0 1 -5 5])
        axis square
        set(gca,'fontSize',16)
    end
end

figure(3)

for i = 1:5
    for j = 1:5
        subplot(5,5,i+5*(j-1))
        
        cma = [nanmean(Cm2,1) nanmean(Cm3,1) nanmean(Cm4,1)];
        
        %plot(nanmean(Cm1,1),nanmean(squeeze(Lambda21(i,j,:,:)),1),'k.')
        hold on
        plot(nanmean(Cm2,1),sqrt(abs(nanmean(squeeze(Lambda22(i,j,:,:)),1))),'g.')
        plot(nanmean(Cm3,1),sqrt(abs(nanmean(squeeze(Lambda23(i,j,:,:)),1))),'b.')
        plot(nanmean(Cm4,1),sqrt(abs(nanmean(squeeze(Lambda24(i,j,:,:)),1))),'r.')
        
        set(gca,'yscale','log');
        axis square
        set(gca,'fontSize',16)
    end
end

% subplot(5,5,1)
% ylabel('0','interpreter','latex')
% subplot(5,5,6)
% ylabel('1','interpreter','latex')
% subplot(5,5,11)
% ylabel('2','interpreter','latex')
% subplot(5,5,16)
% ylabel('3','interpreter','latex')
% subplot(5,5,21)
% ylabel('4','interpreter','latex')
% xlabel('0','interpreter','latex')
% subplot(5,5,22)
% xlabel('1','interpreter','latex')
% subplot(5,5,23)
% xlabel('2','interpreter','latex')
% subplot(5,5,24)
% xlabel('3','interpreter','latex')
% subplot(5,5,25)
% xlabel('4','interpreter','latex')

sgtitle('$M_{i,j}$','interpreter','latex')
% 
% figure(3)
% 
% for i = 1:5
%     for j = 1:5
%         subplot(5,5,i+5*(j-1))
%         
%         cma = [nanmean(Cm2,1) nanmean(Cm3,1) nanmean(Cm4,1)];
%         lm = [nanmean(squeeze(N2(i,j,:,:)),1) nanmean(squeeze(N3(i,j,:,:)),1) nanmean(squeeze(N4(i,j,:,:)),1)];
%         
%         plot(cma,lm,'.','MarkerSize',4)
%         
% %         plot(nanmean(Cm1,1),nanmean(squeeze(N1(i,j,:,:)),1),'+')
% %         hold on
% %         plot(nanmean(Cm2,1),nanmean(squeeze(N2(i,j,:,:)),1),'+')
% %         plot(nanmean(Cm3,1),nanmean(squeeze(N3(i,j,:,:)),1),'+')
% %         plot(nanmean(Cm4,1),nanmean(squeeze(N4(i,j,:,:)),1),'+')
%         
%         axis([0 1 -5 5])
%         axis square
%         set(gca,'fontSize',16)
% 
%     end
% end
% 
% subplot(5,5,1)
% ylabel('0','interpreter','latex')
% subplot(5,5,6)
% ylabel('1','interpreter','latex')
% subplot(5,5,11)
% ylabel('2','interpreter','latex')
% subplot(5,5,16)
% ylabel('3','interpreter','latex')
% subplot(5,5,21)
% ylabel('4','interpreter','latex')
% xlabel('0','interpreter','latex')
% subplot(5,5,22)
% xlabel('1','interpreter','latex')
% subplot(5,5,23)
% xlabel('2','interpreter','latex')
% subplot(5,5,24)
% xlabel('3','interpreter','latex')
% subplot(5,5,25)
% xlabel('4','interpreter','latex')
% 
% sgtitle('$N_{i,j}$','interpreter','latex')