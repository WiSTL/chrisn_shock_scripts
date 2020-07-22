n = 3;

[cm1,um1,vm1,Re1,eps1,scalar_eps1,Rm,Etam] = chris_exp_1D_structure_fn(n,1);
[cm2,um2,vm2,Re2,eps2,scalar_eps2,~,~] = chris_exp_1D_structure_fn(n,0);
[cm3,um3,vm3,Re3,eps3,scalar_eps3,~,~] = chris_exp_1D_structure_fn(n,0);
[cm4,um4,vm4,Re4,eps4,scalar_eps4,~,~] = chris_exp_1D_structure_fn(n,0);

a1 = nanmean(cm1);
a2 = nanmean(cm2);
a3 = nanmean(cm3);
a4 = nanmean(cm4);

b1 = nanmean(um1);
b2 = nanmean(um2);
b3 = nanmean(um3);
b4 = nanmean(um4);

c1 = nanmean(vm1);
c2 = nanmean(vm2);
c3 = nanmean(vm3);
c4 = nanmean(vm4);


figure(1)
plot(Rm(1,:),a1)
hold on
plot(Rm(1,:),a2)
plot(Rm(1,:),a3)
plot(Rm(1,:),a4)

set(gca,'xScale','log')
set(gca,'yScale','log')

figure(2)
plot(Rm(1,:),b1)
hold on
plot(Rm(1,:),b2)
plot(Rm(1,:),b3)
plot(Rm(1,:),b4)

set(gca,'xScale','log')
set(gca,'yScale','log')

figure(3)
plot(Rm(1,:),c1)
hold on
plot(Rm(1,:),c2)
plot(Rm(1,:),c3)
plot(Rm(1,:),c4)

set(gca,'xScale','log')
set(gca,'yScale','log')