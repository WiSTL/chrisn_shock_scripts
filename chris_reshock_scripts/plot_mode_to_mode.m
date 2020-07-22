% [Tuub1,Tccb1,Tuum1,Tccm1,Tuut1,Tcct1] = chris_exp_mode_to_mode(1);
% [Tuub2,Tccb2,Tuum2,Tccm2,Tuut2,Tcct2] = chris_exp_mode_to_mode(0);
% [Tuub3,Tccb3,Tuum3,Tccm3,Tuut3,Tcct3] = chris_exp_mode_to_mode(0);
% [Tuub4,Tccb4,Tuum4,Tccm4,Tuut4,Tcc41] = chris_exp_mode_to_mode(0);

figure
%%%%%%%
% Tuu1
%%%%%%%
subplot(4,3,1)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tuut1,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,2)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tuum1,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,3)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tuub1,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

%%%%%%%
% Tuu2
%%%%%%%
subplot(4,3,4)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tuut2,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,5)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tuum2,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,6)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tuub2,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

%%%%%%%
% Tuu3
%%%%%%%
subplot(4,3,7)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tuut3,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,8)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tuum3,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,9)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tuub3,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

%%%%%%%
% Tuu4
%%%%%%%
subplot(4,3,10)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tuut4,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,11)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tuum4,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,12)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tuub4,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
%%%%%%%
% Tcc1
%%%%%%%
subplot(4,3,1)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tcct1,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,2)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tccm1,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,3)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tccb1,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

%%%%%%%
% Tuu2
%%%%%%%
subplot(4,3,4)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tcct2,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,5)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tccm2,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,6)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tccb2,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

%%%%%%%
% Tuu3
%%%%%%%
subplot(4,3,7)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tcct3,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,8)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tccm3,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,9)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tccb3,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

%%%%%%%
% Tcc4
%%%%%%%
subplot(4,3,10)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tcc41,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,11)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tccm4,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

subplot(4,3,12)
pcolor([0.02:0.02:0.4],[0.02:0.02:0.4],nanmean(Tccb4,3))

axis square
caxis([-0.005 0.005])
shading flat; cmocean('balance')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(4,1,1)
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tcct1,3))/20))
hold on
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tccm1,3))/20))
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tccb1,3))/20))

subplot(4,1,2)
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tcct2,3))/20))
hold on
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tccm2,3))/20))
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tccb2,3))/20))

subplot(4,1,3)
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tcct3,3))/20))
hold on
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tccm3,3))/20))
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tccb3,3))/20))

subplot(4,1,4)
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tcc41,3))/20))
hold on
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tccm4,3))/20))
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tccb4,3))/20))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(4,1,1)
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tuut1,3))/20))
hold on
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tuum1,3))/20))
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tuub1,3))/20))

subplot(4,1,2)
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tuut2,3))/20))
hold on
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tuum2,3))/20))
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tuub2,3))/20))

subplot(4,1,3)
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tuut3,3))/20))
hold on
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tuum3,3))/20))
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tuub3,3))/20))

subplot(4,1,4)
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tuut4,3))/20))
hold on
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tuum4,3))/20))
plot([0.02:0.02:0.4],cumtrapz(trapz(nanmean(Tuub4,3))/20))

