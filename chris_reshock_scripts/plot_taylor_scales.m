
[d1,tau1] = chris_exp_Taylor_Scales(1,[-2:0.001:2]);
[d2,tau2] = chris_exp_Taylor_Scales(0,[-1.6:0.001:1]);
[d3,tau3] = chris_exp_Taylor_Scales(0,[-1.2:0.001:0.8]);
[d4,tau4] = chris_exp_Taylor_Scales(0,[-0.8:0.001:0.5]);

%%%%%%%%%%%%%%%%%%%%%%

figure(1)
n=1;
m=2;
l=12;

plot(d1(n,:)./d1(l,:),d1(m,:)./d1(l,:),'k+')
hold on
plot(d2(n,:)./d2(l,:),d2(m,:)./d2(l,:),'bs')
plot(d3(n,:)./d3(l,:),d3(m,:)./d3(l,:),'gx')
plot(d4(n,:)./d4(l,:),d4(m,:)./d4(l,:),'rv')

set(gca,'xscale','log')
set(gca,'yscale','log')
axis square
axis([0.005 0.1 0.005 0.1])


figure(2)
n=3;
m=4;
l=12;

plot(d1(n,:)./d1(l,:),d1(m,:)./d1(l,:),'k+')
hold on
plot(d2(n,:)./d2(l,:),d2(m,:)./d2(l,:),'bs')
plot(d3(n,:)./d3(l,:),d3(m,:)./d3(l,:),'gx')
plot(d4(n,:)./d4(l,:),d4(m,:)./d4(l,:),'rv')

figure(3)
n=5;
m=6;
l=13;

plot(d1(n,:)./d1(l,:),d1(m,:)./d1(l,:),'k+')
hold on
plot(d2(n,:)./d2(l,:),d2(m,:)./d2(l,:),'bs')
plot(d3(n,:)./d3(l,:),d3(m,:)./d3(l,:),'gx')
plot(d4(n,:)./d4(l,:),d4(m,:)./d4(l,:),'rv')

figure(4)
n=7;
m=8;

plot(d1(n,:)./sqrt(d1(12,:).*d1(13,:)),d1(m,:)./sqrt(d1(12,:).*d1(13,:)),'k+')
hold on
plot(d2(n,:)./sqrt(d2(12,:).*d2(13,:)),d2(m,:)./sqrt(d2(12,:).*d2(13,:)),'bs')
plot(d3(n,:)./sqrt(d3(12,:).*d3(13,:)),d3(m,:)./sqrt(d3(12,:).*d3(13,:)),'gx')
plot(d4(n,:)./sqrt(d4(12,:).*d4(13,:)),d4(m,:)./sqrt(d4(12,:).*d4(13,:)),'rv')

figure(5)
n=9;
m=10;

plot(d1(n,:)./sqrt(d1(12,:).*d1(13,:)),d1(m,:)./sqrt(d1(12,:).*d1(13,:)),'k+')
hold on
plot(d2(n,:)./sqrt(d2(12,:).*d2(13,:)),d2(m,:)./sqrt(d2(12,:).*d2(13,:)),'bs')
plot(d3(n,:)./sqrt(d3(12,:).*d3(13,:)),d3(m,:)./sqrt(d3(12,:).*d3(13,:)),'gx')
plot(d4(n,:)./sqrt(d4(12,:).*d4(13,:)),d4(m,:)./sqrt(d4(12,:).*d4(13,:)),'rv')

figure(6)
n=12;
m=13;

plot(d1(n,:),d1(m,:),'k+')
hold on
plot(d2(n,:),d2(m,:),'bs')
plot(d3(n,:),d3(m,:),'gx')
plot(d4(n,:),d4(m,:),'rv')

figure(7)
n=12;
m=13;

plot(tau1,d1(m,:)./d1(n,:),'k+')
hold on
plot(tau2,d2(m,:)./d2(n,:),'bs')
plot(tau3,d3(m,:)./d3(n,:),'gx')
plot(tau4,d4(m,:)./d4(n,:),'rv')

figure(8)
n=6;
m=5;

plot(tau1,d1(m,:)./d1(n,:),'k+')
hold on
plot(tau2,d2(m,:)./d2(n,:),'bs')
plot(tau3,d3(m,:)./d3(n,:),'gx')
plot(tau4,d4(m,:)./d4(n,:),'rv')

figure(9)
n=11;

plot(tau1,d1(n,:),'k+')
hold on
plot(tau2,d2(n,:),'bs')
plot(tau3,d3(n,:),'gx')
plot(tau4,d4(n,:),'rv')


figure(10)
n=13;
m=15;
l=12;
p=14;

plot(tau1,(d1(m,:)./d1(n,:)).*(d1(l,:)./d1(p,:)),'k+')
hold on
plot(tau2,(d2(m,:)./d2(n,:)).*(d2(l,:)./d2(p,:)),'bs')
plot(tau3,(d3(m,:)./d3(n,:)).*(d3(l,:)./d3(p,:)),'gx')
plot(tau4,(d4(m,:)./d4(n,:)).*(d4(l,:)./d4(p,:)),'rv')

figure(11)
m=13;
n=15;
p=12;
l=14;

plot(sqrt(3200/3)*(d1(m,:)./d1(n,:)),sqrt(3200/3)*(d1(p,:)./d1(l,:)),'k+')
hold on
plot(sqrt(3200/3)*(d2(m,:)./d2(n,:)),sqrt(3200/3)*(d2(p,:)./d2(l,:)),'bs')
plot(sqrt(3200/3)*(d3(m,:)./d3(n,:)),sqrt(3200/3)*(d3(p,:)./d3(l,:)),'gx')
plot(sqrt(3200/3)*(d4(m,:)./d4(n,:)),sqrt(3200/3)*(d4(p,:)./d4(l,:)),'rv')

figure(12)
n=5;
m=13;
l=6;
p=13;

plot(sqrt(3200/3)*(d1(m,:)./d1(n,:)),sqrt(3200/3)*(d1(p,:)./d1(l,:)),'k+')
hold on
plot(sqrt(3200/3)*(d2(m,:)./d2(n,:)),sqrt(3200/3)*(d2(p,:)./d2(l,:)),'bs')
plot(sqrt(3200/3)*(d3(m,:)./d3(n,:)),sqrt(3200/3)*(d3(p,:)./d3(l,:)),'gx')
plot(sqrt(3200/3)*(d4(m,:)./d4(n,:)),sqrt(3200/3)*(d4(p,:)./d4(l,:)),'rv')

figure(13)
n=1;
m=4;
l=12;

plot(d1(n,:)./d1(l,:),d1(m,:)./d1(l,:),'k+')
hold on
plot(d2(n,:)./d2(l,:),d2(m,:)./d2(l,:),'bs')
plot(d3(n,:)./d3(l,:),d3(m,:)./d3(l,:),'gx')
plot(d4(n,:)./d4(l,:),d4(m,:)./d4(l,:),'rv')