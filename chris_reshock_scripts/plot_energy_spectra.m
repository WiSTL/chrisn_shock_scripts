[cm1,em1,un1,vn1,Km,Etam] = chris_exp_1D_spectra_spatial_2018a(1);
[cm2,em2,un2,vn2,~,~] = chris_exp_1D_spectra_spatial_2018a(0);
[cm3,em3,un3,vn3,~,~] = chris_exp_1D_spectra_spatial_2018a(0);
[cm4,em4,un4,vn4,~,~] = chris_exp_1D_spectra_spatial_2018a(0);

a1 = nanmean(abs(cm1));
a2 = nanmean(abs(cm2));
a3 = nanmean(abs(cm3));
a4 = nanmean(abs(cm4));

b1 = nanmean(abs(em1));
b2 = nanmean(abs(em2));
b3 = nanmean(abs(em3));
b4 = nanmean(abs(em4));

k = Km(1,:);
figure(5)
plot(Km(1,:),(k.^(2)).*a1./max(a1))
hold on
plot(Km(1,:),(k.^(2)).*a2./max(a2))
plot(Km(1,:),(k.^(2)).*a3./max(a3))
plot(Km(1,:),(k.^(2)).*a4./max(a4))

set(gca,'xScale','log')
set(gca,'yScale','log')

figure(6)
plot(Km(1,:),(k.^(1.127)).*b1./max(b1))
hold on
plot(Km(1,:),(k.^(1.127)).*b2./max(b2))
plot(Km(1,:),(k.^(1.127)).*b3./max(b3))
plot(Km(1,:),(k.^(1.127)).*b4./max(b4))

set(gca,'xScale','log')
set(gca,'yScale','log')
