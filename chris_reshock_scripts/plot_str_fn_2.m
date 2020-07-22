% [zeta_c_n1,zeta_u_n1,zeta_w_n1,Sc21,Su21,Sw21,rx1,Scn1,Sun1,Swn1] = chris_exp_str_fn(1);
% [zeta_c_n2,zeta_u_n2,zeta_w_n2,Sc22,Su22,Sw22,rx2,Scn2,Sun2,Swn2] = chris_exp_str_fn(0);
% [zeta_c_n3,zeta_u_n3,zeta_w_n3,Sc23,Su23,Sw23,rx3,Scn3,Sun3,Swn3] = chris_exp_str_fn(0);
% [zeta_c_n4,zeta_u_n4,zeta_w_n4,Sc24,Su24,Sw24,rx4,Scn4,Sun4,Swn4] = chris_exp_str_fn(0);


figure(1)

errorbar([2:10],nanmean(imgaussfilt(zeta_c_n1(2:10,:),0.9),2),nanvar(imgaussfilt(zeta_c_n1(2:10,:),0.9)'),nanvar(imgaussfilt(zeta_c_n1(2:10,:),0.9)'),'k+')
hold on
errorbar([2:10],nanmean(imgaussfilt(zeta_c_n2(2:10,:),0.9),2),nanvar(imgaussfilt(zeta_c_n2(2:10,:),0.9)'),nanvar(imgaussfilt(zeta_c_n2(2:10,:),0.9)'),'bs')
errorbar([2:10],nanmean(imgaussfilt(zeta_c_n3(2:10,:),0.9),2),nanvar(imgaussfilt(zeta_c_n3(2:10,:),0.9)'),nanvar(imgaussfilt(zeta_c_n3(2:10,:),0.9)'),'gx')
errorbar([2:10],nanmean(imgaussfilt(zeta_c_n4(2:10,:),0.9),2),nanvar(imgaussfilt(zeta_c_n4(2:10,:),0.9)'),nanvar(imgaussfilt(zeta_c_n4(2:10,:),0.9)'),'rv')


figure(2)

errorbar([2:10],nanmean(imgaussfilt(zeta_u_n1(2:10,:),0.9),2),nanvar(imgaussfilt(zeta_u_n1(2:10,:),0.9)'),nanvar(imgaussfilt(zeta_u_n1(2:10,:),0.9)'),'k+')
hold on
errorbar([2:10],nanmean(imgaussfilt(zeta_u_n2(2:10,:),0.9),2),nanvar(imgaussfilt(zeta_u_n2(2:10,:),0.9)'),nanvar(imgaussfilt(zeta_u_n2(2:10,:),0.9)'),'bs')
errorbar([2:10],nanmean(imgaussfilt(zeta_u_n3(2:10,:),0.9),2),nanvar(imgaussfilt(zeta_u_n3(2:10,:),0.9)'),nanvar(imgaussfilt(zeta_u_n3(2:10,:),0.9)'),'gx')
errorbar([2:10],nanmean(imgaussfilt(zeta_u_n4(2:10,:),0.9),2),nanvar(imgaussfilt(zeta_u_n4(2:10,:),0.9)'),nanvar(imgaussfilt(zeta_u_n4(2:10,:),0.9)'),'rv')


figure(3)

errorbar([2:10],nanmean(imgaussfilt(zeta_w_n1(2:10,:),0.9),2),nanvar(imgaussfilt(zeta_w_n1(2:10,:),0.9)'),nanvar(imgaussfilt(zeta_w_n1(2:10,:),0.9)'),'k+')
hold on
errorbar([2:10],nanmean(imgaussfilt(zeta_w_n2(2:10,:),0.9),2),nanvar(imgaussfilt(zeta_w_n2(2:10,:),0.9)'),nanvar(imgaussfilt(zeta_w_n2(2:10,:),0.9)'),'bs')
errorbar([2:10],nanmean(imgaussfilt(zeta_w_n3(2:10,:),0.9),2),nanvar(imgaussfilt(zeta_w_n3(2:10,:),0.9)'),nanvar(imgaussfilt(zeta_w_n3(2:10,:),0.9)'),'gx')
errorbar([2:10],nanmean(imgaussfilt(zeta_w_n4(2:10,:),0.9),2),nanvar(imgaussfilt(zeta_w_n4(2:10,:),0.9)'),nanvar(imgaussfilt(zeta_w_n4(2:10,:),0.9)'),'rv')




% figure(1)
% chris_plot_mean_std([2:10],(zeta_c_n1(2:10,:)),'k',10);
% hold on
% chris_plot_mean_std([2:10],(zeta_c_n2(2:10,:)),'b',10);
% hold on
% chris_plot_mean_std([2:10],(zeta_c_n3(2:10,:)),'g',10);
% hold on
% chris_plot_mean_std([2:10],(zeta_c_n4(2:10,:)),'r',10);
% 
% figure(2)
% chris_plot_mean_std([2:10],(zeta_u_n1(2:10,:)),'k',10);
% hold on
% chris_plot_mean_std([2:10],(zeta_u_n2(2:10,:)),'b',10);
% hold on
% chris_plot_mean_std([2:10],(zeta_u_n3(2:10,:)),'g',10);
% hold on
% chris_plot_mean_std([2:10],(zeta_u_n4(2:10,:)),'r',10);
% 
% figure(3)
% chris_plot_mean_std([2:10],(zeta_w_n1(2:10,:)),'k',10);
% hold on
% chris_plot_mean_std([2:10],(zeta_w_n2(2:10,:)),'b',10);
% hold on
% chris_plot_mean_std([2:10],(zeta_w_n3(2:10,:)),'g',10);
% hold on
% chris_plot_mean_std([2:10],(zeta_w_n4(2:10,:)),'r',10);