[Pn1,Tn1,Dn1,en1,ddn1,termn1,eta_m1] = chris_exp_equation_terms(1,[-3:0.001:2]);
[Pn2,Tn2,Dn2,en2,ddn2,termn2,eta_m2] = chris_exp_equation_terms(0,[-1.6:0.001:1.6]);
[Pn3,Tn3,Dn3,en3,ddn3,termn3,eta_m3] = chris_exp_equation_terms(0,[-1.2:0.001:1.2]);

figure(1)
chris_plot_mean_std(eta_m1,squeeze(Pn1(1,:,:))','k');
hold on
chris_plot_mean_std(eta_m2,squeeze(Pn2(1,:,:))','b');
hold on
chris_plot_mean_std(eta_m3,squeeze(Pn3(1,:,:))','g');

figure(2)
chris_plot_mean_std(eta_m1,squeeze(Tn1(1,:,:))','k');
hold on
chris_plot_mean_std(eta_m2,squeeze(Tn2(1,:,:))','b');
hold on
chris_plot_mean_std(eta_m3,squeeze(Tn3(1,:,:))','g');

figure(3)
chris_plot_mean_std(eta_m1,squeeze(Pn1(2,:,:))','k');
hold on
chris_plot_mean_std(eta_m2,squeeze(Pn2(2,:,:))','b');
hold on
chris_plot_mean_std(eta_m3,squeeze(Pn3(2,:,:))','g');

figure(4)
chris_plot_mean_std(eta_m1,squeeze(Tn1(2,:,:))','k');
hold on
chris_plot_mean_std(eta_m2,squeeze(Tn2(2,:,:))','b');
hold on
chris_plot_mean_std(eta_m3,squeeze(Tn3(2,:,:))','g');

figure(5)
chris_plot_mean_std(eta_m1,squeeze(Pn1(3,:,:))','k');
hold on
chris_plot_mean_std(eta_m2,squeeze(Pn2(3,:,:))','b');
hold on
chris_plot_mean_std(eta_m3,squeeze(Pn3(3,:,:))','g');

figure(6)
chris_plot_mean_std(eta_m1,squeeze(Tn1(3,:,:))','k');
hold on
chris_plot_mean_std(eta_m2,squeeze(Tn2(3,:,:))','b');
hold on
chris_plot_mean_std(eta_m3,squeeze(Tn3(3,:,:))','g');

figure(7)
chris_plot_mean_std(eta_m1,squeeze(Pn1(4,:,:))','k');
hold on
chris_plot_mean_std(eta_m2,squeeze(Pn2(4,:,:))','b');
hold on
chris_plot_mean_std(eta_m3,squeeze(Pn3(4,:,:))','g');

figure(8)
chris_plot_mean_std(eta_m1,squeeze(Tn1(4,:,:))','k');
hold on
chris_plot_mean_std(eta_m2,squeeze(Tn2(4,:,:))','b');
hold on
chris_plot_mean_std(eta_m3,squeeze(Tn3(4,:,:))','g');

figure(9)
chris_plot_mean_std(eta_m1,squeeze(Pn1(5,:,:))','k');
hold on
chris_plot_mean_std(eta_m2,squeeze(Pn2(5,:,:))','b');
hold on
chris_plot_mean_std(eta_m3,squeeze(Pn3(5,:,:))','g');

figure(10)
chris_plot_mean_std(eta_m1,squeeze(Tn1(5,:,:))','k');
hold on
chris_plot_mean_std(eta_m2,squeeze(Tn2(5,:,:))','b');
hold on
chris_plot_mean_std(eta_m3,squeeze(Tn3(5,:,:))','g');