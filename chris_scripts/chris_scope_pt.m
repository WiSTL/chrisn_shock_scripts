function [xloc,pt,pt_t,t,xt,ts,trs,Ws,Wt,Wr,Ms,Wsv] = chris_scope_pt(fn,vars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Reads lvm files from the scope program,
%   normalises pressure transducer data
%   smooths data and calculates derivatives and
%   estimates wave speeds.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pt=dlmread(fn,'\t',[23 1 40022 12]);

for i = 1:12
pt(:,i) = pt(:,i)/max(smooth(pt(:,i),50));
end

% pt2 = pt(:,10);
% pt(:,10) = pt(:,9);
% pt(:,9) = pt2;

pt(1:100,:) = 0;
pt = pt(:,[1 2 3 4 5 6 7 12 9 10 11 8]);

xloc = [0.501435 1.7145 2.9464 3.697 4.4021 4.628 4.628 4.929188 6.9596125-1.31445 6.9596125-0.75565 6.9596125-0.24765 6.9596125]; %6.7072

ltot = xloc(end);
window2_top = ltot - 0.4064;
window2_bottom = ltot - 0.1778;

window3_top = ltot - 0.620713;
window3_bottom = ltot - 0.392113;

window4_top = ltot - 1.0414;
window4_bottom = ltot - 0.9398;

window1_top = ltot - 2.0955;
window1_bottom = ltot - 1.9939;

%%% old measurments
%%% xloc = [0.5017 1.715 2.927 3.697 4.3066 4.5225 4.5225 4.8273 5.221 5.5258 6.6053 6.8339];

dt = 1*10^(-6);

t = dt:dt:40000*dt;



xt = pt*0;

for i = 1:12
    pt_t(:,i) = smooth(chris_derivative(smooth(pt(:,i),100),1*10^(-6)),100);
    
    a = pt_t(:,i);
    
    pt_t(:,i) = pt_t(:,i)/max(a(:));
    
    xt(:,i)=a>0.15*max(a);
    
    tss = find(xt(:,i)==1);
    ts(i) = t(tss(1));
    
    a = t(tss);
    
    trss = find(a>0.35*max(t));
    if numel(trss)>0
        trs(i) = a(trss(1));
    else
        trs(i) = nan;
    end
    
end

ts(5) = vars(14);
ts(6) = vars(13);
ts(7) = vars(12);
ts(8) = vars(11);
ts(9) = vars(15);
ts(10) = vars(10);
ts(11) = vars(7);
trs(11) = vars(8);
ts(12) = vars(9);
trs(12) = vars(9);

%%%
% Incident
%%%

x_i = xloc(6:9);
t_i = ts(6:9);

% figure(5)
% plot(t_i,x_i)

p = polyfit(t_i,x_i,1);

Ws = p(1);


%%%
% transmitted
%%%

x1 = 0.3302;
x2 = 0.2286;

dt = ts(10)-ts(9);

Wt1 = x2/(dt - (x1/Ws))

x_t = xloc([10:12]);
t_t = ts([10:12]);

ptr = polyfit(t_t,x_t,1);

Wt = ptr(1);

%%%
% reflected
%%%

x_r = [xloc(11),xloc(12)];
t_r = [trs(11),trs(12)];

pr = polyfit(t_r,x_r,1);

Wr = pr(1);


Wsv(1) = nan;

for i = 2:12
    if i==7
        Wsv(i) = nan;
    else
        Wsv(i) = 1/nanmean(diff(ts(i-1:i))./diff(xloc(i-1:i)));
    end
end

% Ws = (xloc(10)-xloc(4))/(ts(10)-ts(4));


%%%
% EES values
%%%

% a1 = 1011;%He
% a1 = 676.4; %He + 7.5% acetone; 
%a1 = 377.3; %He+12.18%
%a1 = 376.5; %He + 33.32% acetone; 

a1 = vars(1);
ua3 = vars(2);
Wc2 = vars(3);
ua4 = vars(4);

triggerpt = vars(16);


%%%
% End EES values
%%%

%%%
% Interface Location
%%%

triggerpt_delay = vars(5);

adj = vars(6);

x_int_start = 5.97;

t_int_start = ts(8) + (x_int_start-xloc(8))/Ws;
t_int_final = 0.01+triggerpt_delay;

x_int_final = x_int_start+ua3*(t_int_final-adj - t_int_start);

x_int = [x_int_start,x_int_start,x_int_final];
t_int = [0,t_int_start,t_int_final-adj];

%%%
% reshock Location
%%%

x_rs_start = ltot;

t_rs_start = trs(12);
t_rs_final = t_int_final;

x_rs_final = x_rs_start+Wc2*(t_rs_final - t_rs_start);

x_rs = [x_rs_start,x_rs_final];
t_rs = [t_rs_start,t_rs_final];

%%%
% reshocked interface Location
%%%

x_rsi_start = x_int_final;

t_rsi_start = t_int_final-adj;
t_rsi_final = t_int_final;

x_rsi_final = x_rsi_start+ua4*(t_rsi_final - t_rsi_start);

x_rsi = [x_rsi_start,x_rsi_final];
t_rsi = [t_rsi_start,t_rsi_final];

% ts = ts+0.0015;
% trs = trs+0.0015;
% t_rs = t_rs+0.0015; %Run12 timing


g = 1.67;

Ms = Ws/a1;

up = (2*a1/(g+1))*(Ms-(1/Ms));

tcs = t(1:end/7);
xcs = up*tcs;

trs(6:7) = nan;

% dt = -(0.0097-0.0088);
% figure(1)
% hold on
% plot(xloc,ts+dt,'+')
% hold on
% plot(xloc,trs+dt,'+')
% plot(xloc,ts+dt)
% plot(xloc,trs+dt)
% plot(6.65,t_int_final+dt,'rO')
% plot(xcs,tcs+0.003)
% plot(x_int,t_int)
% plot(x_rs,t_rs+dt)
% plot(x_rsi,t_rsi)
% grid on
% xlabel('x [m]')
% ylabel('t [s]')
% legend('incident shock path','reflected shock path')
% title('shock trajectories')

t_int*1000
t_int_final*1000

figure(1)
hold on
plot(xloc,ts*1000,'+')
hold on
plot(xloc,trs*1000,'+')
plot(xloc,ts*1000)
plot(xloc,trs*1000)
plot(6.65,t_int_final*1000,'rO')
plot(xcs,(tcs+0.003)*1000)
plot(x_int,t_int*1000)
plot(x_rs,t_rs*1000)
plot(x_rsi,t_rsi*1000)
grid on
xlabel('x [m]')
ylabel('t [ms]')
legend('incident shock path','reflected shock path')
title('shock trajectories')

win2 = rectangle('position',[window2_top 0 abs(window2_top-window2_bottom)  max(max(t*1000))],'EdgeColor','k','LineStyle','--');
win3 = rectangle('position',[window3_top 0 abs(window3_top-window3_bottom)  max(max(t*1000))],'EdgeColor','k','LineStyle','--');
win4 = rectangle('position',[window4_top 0 abs(window4_top-window4_bottom)  max(max(t*1000))],'EdgeColor','k','LineStyle','--');

[X,T]=meshgrid(xloc,t);
figure(2)
hold off
pcolor(X,T*1000,pt), shading flat
hold on
plot(xloc,ts*1000,'+')
hold on
plot(xloc,trs*1000,'+')
plot(xloc,ts*1000)
plot(xloc,trs*1000)
plot(6.65,t_int_final*1000,'rO')
grid on
title('Pressure X-T diagram')
xlabel('x [m]')
ylabel('t [ms]')

xloc2 = xloc;
pt2 = pt;

xloc2(9)=[];
xloc2(9)=[];
pt2(:,9)=[];
pt2(:,9)=[];

pt2(:,7)=[];
xloc2(7)=[];

figure(3)
pcolor(t*1000,xloc2,pt2')
shading interp
hold on
plot(t_int*1000,x_int,'k','LineWidth',2)
plot(t_rsi*1000,x_rsi,'k','LineWidth',2)
hold off
view([90 -90])
title('Pressure X-T diagram')
xlabel('t [ms]')
ylabel('x [m]')

xloc3 = xloc;
pt3 = pt_t;

xloc3(9)=[];
xloc3(9)=[];
pt3(:,9)=[];
pt3(:,9)=[];

pt3(:,7)=[];
xloc3(7)=[];

figure(4)
pcolor(t*1000,xloc3,pt3')
shading interp
hold on
plot(t_int*1000,x_int,'k','LineWidth',2)
plot(t_rsi*1000,x_rsi,'k','LineWidth',2)
hold off
view([90 -90])
title('Pressure X-T diagram')
xlabel('t [ms]')
ylabel('x [m]')


end

