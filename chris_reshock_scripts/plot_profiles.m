
[up2n1,upvpn1,apupn1,apvpn1,vp2n1,vp3n1,upupvpn1,upvpdvpn1,vpvpdupn1,cn1,cpn1,un1,vn1,xin1,udcn1,vdcn1,t11n1,t12n1,t22n1,bn1,t11Rn1,t12Rn1,t22Rn1,rvn1,run1,eresn1,esgsn1,b11Rn1,b12Rn1,b22Rn1,q2n1,eta_m1,Re1] = chris_exp_profiles(1,[-2:0.001:2]);
[up2n2,upvpn2,apupn2,apvpn2,vp2n2,vp3n2,upupvpn2,upvpdvpn2,vpvpdupn2,cn2,cpn2,un2,vn2,xin2,udcn2,vdcn2,t11n2,t12n2,t22n2,bn2,t11Rn2,t12Rn2,t22Rn2,rvn2,run2,eresn2,esgsn2,b11Rn2,b12Rn2,b22Rn2,q2n2,eta_m2,Re2] = chris_exp_profiles(0,[-1.6:0.001:1]);
[up2n3,upvpn3,apupn3,apvpn3,vp2n3,vp3n3,upupvpn3,upvpdvpn3,vpvpdupn3,cn3,cpn3,un3,vn3,xin3,udcn3,vdcn3,t11n3,t12n3,t22n3,bn3,t11Rn3,t12Rn3,t22Rn3,rvn3,run3,eresn3,esgsn3,b11Rn3,b12Rn3,b22Rn3,q2n3,eta_m3,Re3] = chris_exp_profiles(0,[-1.2:0.001:0.8]);
[up2n4,upvpn4,apupn4,apvpn4,vp2n4,vp3n4,upupvpn4,upvpdvpn4,vpvpdupn4,cn4,cpn4,un4,vn4,xin4,udcn4,vdcn4,t11n4,t12n4,t22n4,bn4,t11Rn4,t12Rn4,t22Rn4,rvn4,run4,eresn4,esgsn4,b11Rn4,b12Rn4,b22Rn4,q2n4,eta_m4,Re4] = chris_exp_profiles(0,[-0.8:0.001:0.5]);

figure(1)
chris_plot_mean_std(eta_m1,(cn1)','k',2);
hold on
chris_plot_mean_std(eta_m2,(cn2)','b',2);
hold on
chris_plot_mean_std(eta_m3,(cn3)','g',2);
hold on
chris_plot_mean_std(eta_m4,(cn4)','r',2);

figure(2)
chris_plot_mean_std(eta_m1,(cpn1)','k',2);
hold on
chris_plot_mean_std(eta_m2,(cpn2)','b',2);
hold on
chris_plot_mean_std(eta_m3,(cpn3)','g',2);
hold on
chris_plot_mean_std(eta_m4,(cpn4)','r',2);

figure(3)
chris_plot_mean_std(eta_m1,(vn1)','k',2);
hold on
chris_plot_mean_std(eta_m2,(vn2)','b',2);
hold on
chris_plot_mean_std(eta_m3,(vn3)','g',2);
hold on
chris_plot_mean_std(eta_m4,(vn4)','r',2);


figure(4)
chris_plot_mean_std(eta_m1,(un1)','k',2);
hold on
chris_plot_mean_std(eta_m2,(un2)','b',2);
hold on
chris_plot_mean_std(eta_m3,(un3)','g',2);
hold on
chris_plot_mean_std(eta_m4,(un4)','r',2);

figure(5)
chris_plot_mean_std(eta_m1,(vp2n1)','k',2);
hold on
chris_plot_mean_std(eta_m2,(vp2n2)','b',2);
hold on
chris_plot_mean_std(eta_m3,(vp2n3)','g',2);
hold on
chris_plot_mean_std(eta_m4,(vp2n4)','r',2);


figure(6)
chris_plot_mean_std(eta_m1,(vp3n1)','k',2);
hold on
chris_plot_mean_std(eta_m2,(vp3n2)','b',2);
hold on
chris_plot_mean_std(eta_m3,(vp3n3)','g',2);
hold on
chris_plot_mean_std(eta_m4,(vp3n4)','r',2);

figure(7)
chris_plot_mean_std(eta_m1,(b11Rn1)','k',2);
hold on
chris_plot_mean_std(eta_m2,(b11Rn2)','b',2);
hold on
chris_plot_mean_std(eta_m3,(b11Rn3)','g',2);
hold on
chris_plot_mean_std(eta_m4,(b11Rn4)','r',2);

figure(8)
chris_plot_mean_std(eta_m1,(b12Rn1)','k',2);
hold on
chris_plot_mean_std(eta_m2,(b12Rn2)','b',2);
hold on
chris_plot_mean_std(eta_m3,(b12Rn3)','g',2);
hold on
chris_plot_mean_std(eta_m4,(b12Rn4)','r',2);

figure(9)
chris_plot_mean_std(eta_m1,(b22Rn1)','k',2);
hold on
chris_plot_mean_std(eta_m2,(b22Rn2)','b',2);
hold on
chris_plot_mean_std(eta_m3,(b22Rn3)','g',2);
hold on
chris_plot_mean_std(eta_m4,(b22Rn4)','r',2);

figure(10)
chris_plot_mean_std(eta_m1,(rpvpn1)','k',2);
hold on
chris_plot_mean_std(eta_m2,(rpvpn2)','b',2);
hold on
chris_plot_mean_std(eta_m3,(rpvpn3)','g',2);
hold on
chris_plot_mean_std(eta_m4,(rpvpn4)','r',2);

figure(11)
chris_plot_mean_std(eta_m1,(apvpn1)','k',2);
hold on
chris_plot_mean_std(eta_m2,(apvpn2)','b',2);
hold on
chris_plot_mean_std(eta_m3,(apvpn3)','g',2);
hold on
chris_plot_mean_std(eta_m4,(apvpn4)','r',2);