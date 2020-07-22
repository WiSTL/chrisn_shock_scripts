function [Mt,Mr,Mwr,Pt,Pr,Pwr,Ps] = chris_shocked_interface_He_Ar(Ms)

rho_he = 0.164;
rho_ar = 1.449;

T_atm = 300;

Mr = 1.377 - 0.9752*exp(-0.9558*Ms);
Mt = 1.529*Ms - 0.5572;
Mwr = 2.15 - 3.846*exp(-1.207*Ms);

[Ps,rhos,Ts,P0s,T0s,M2s] = moving_shock(Ms,1.503);
[Pr,rhor,Tr,P0r,T0r,M2r] = moving_shock(Mr,1.503);
[Pt,rhot,Tt,P0t,T0t,M2t] = moving_shock(Mt,1.664);
[Pwr,rhowr,Twr,P0wr,T0wr,M2wr] = moving_shock(Mwr,1.664);

figure
plot(Ms,rhos.*rhor*rho_he)
hold on
plot(Ms,Ts.*Tr)
legend('Density','Temperature')
title('region 1')

figure
plot(Ms,rhot*rho_ar)
hold on
plot(Ms,Tt)
legend('Density','Temperature')
title('region 2')





end

