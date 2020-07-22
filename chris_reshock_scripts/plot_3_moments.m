
[un1,up2n1,up3n1,vn1,vp2n1,vp3n1,cn1,cp2n1,cp3n1,upvp2n1,luxn1,luyn1,lvxn1,lvyn1,lcxn1,lcyn1,epsxn1,epsyn1,Rt1,R1,del1,beta1,tau1,eta_m1] = chris_exp_3_moments(1,[-2:0.001:2]);
[un2,up2n2,up3n2,vn2,vp2n2,vp3n2,cn2,cp2n2,cp3n2,upvp2n2,luxn2,luyn2,lvxn2,lvyn2,lcxn2,lcyn2,epsxn2,epsyn2,Rt2,R2,del2,beta2,tau2,eta_m2] = chris_exp_3_moments(0,[-1:0.001:1.6]);
[un3,up2n3,up3n3,vn3,vp2n3,vp3n3,cn3,cp2n3,cp3n3,upvp2n3,luxn3,luyn3,lvxn3,lvyn3,lcxn3,lcyn3,epsxn3,epsyn3,Rt3,R3,del3,beta3,tau3,eta_m3] = chris_exp_3_moments(0,[-0.8:0.001:1.2]);
[un4,up2n4,up3n4,vn4,vp2n4,vp3n4,cn4,cp2n4,cp3n4,upvp2n4,luxn4,luyn4,lvxn4,lvyn4,lcxn4,lcyn4,epsxn4,epsyn4,Rt4,R4,del4,beta4,tau4,eta_m4] = chris_exp_3_moments(0,[-0.5:0.001:0.8]);

figure(1)
chris_plot_mean_std(eta_m1,(cn1)','k');
hold on
chris_plot_mean_std(eta_m2,(cn2)','b');
hold on
chris_plot_mean_std(eta_m3,(cn3)','g');
hold on
chris_plot_mean_std(eta_m4,(cn4)','r');

figure(2)
chris_plot_mean_std(eta_m1,(cp2n1)','k');
hold on
chris_plot_mean_std(eta_m2,(cp2n2)','b');
hold on
chris_plot_mean_std(eta_m3,(cp2n3)','g');
hold on
chris_plot_mean_std(eta_m4,(cp2n4)','r');

figure(3)
chris_plot_mean_std(eta_m1,(cp3n1)','k');
hold on
chris_plot_mean_std(eta_m2,(cp3n2)','b');
hold on
chris_plot_mean_std(eta_m3,(cp3n3)','g');
hold on
chris_plot_mean_std(eta_m4,(cp3n4)','r');

figure(4)
chris_plot_mean_std(eta_m1,(vn1)','k');
hold on
chris_plot_mean_std(eta_m2,(vn2)','b');
hold on
chris_plot_mean_std(eta_m3,(vn3)','g');
hold on
chris_plot_mean_std(eta_m4,(vn4)','r');

figure(5)
chris_plot_mean_std(eta_m1,(vp2n1)','k');
hold on
chris_plot_mean_std(eta_m2,(vp2n2)','b');
hold on
chris_plot_mean_std(eta_m3,(vp2n3)','g');
hold on
chris_plot_mean_std(eta_m4,(vp2n4)','r');


figure(6)
chris_plot_mean_std(eta_m1,(vp3n1)','k');
hold on
chris_plot_mean_std(eta_m2,(vp3n2)','b');
hold on
chris_plot_mean_std(eta_m3,(vp3n3)','g');
hold on
chris_plot_mean_std(eta_m4,(vp3n4)','r');

figure(7)
chris_plot_mean_std(eta_m1,(un1)','k');
hold on
chris_plot_mean_std(eta_m2,(un2)','b');
hold on
chris_plot_mean_std(eta_m3,(un3)','g');
hold on
chris_plot_mean_std(eta_m4,(un4)','r');

figure(8)
chris_plot_mean_std(eta_m1,(up2n1)','k');
hold on
chris_plot_mean_std(eta_m2,(up2n2)','b');
hold on
chris_plot_mean_std(eta_m3,(up2n3)','g');
hold on
chris_plot_mean_std(eta_m4,(up2n4)','r');

figure(9)
chris_plot_mean_std(eta_m1,(up3n1)','k');
hold on
chris_plot_mean_std(eta_m2,(up3n2)','b');
hold on
chris_plot_mean_std(eta_m3,(up3n3)','g');
hold on
chris_plot_mean_std(eta_m4,(up3n4)','r');

figure(10)
chris_plot_mean_std(eta_m1,(luxn1)','k');
hold on
chris_plot_mean_std(eta_m2,(luxn2)','b');
hold on
chris_plot_mean_std(eta_m3,(luxn3)','g');
hold on
chris_plot_mean_std(eta_m4,(luxn4)','r');