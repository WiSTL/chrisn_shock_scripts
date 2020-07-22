function [l_ux2,l_uy2,l_vx2,l_vy2,l_cx2,l_cy2,l_cvx2,l_cvy2,l_cux2,l_cuy2,l_K2,l_C2] = chris_taylor_microscale(u,v,C,dpxdx)

u = u-mean2(u);
v = v-mean2(v);

% [Rx,rx,Ry,ry] = chris_autocorrelation(u,u);
% 
% ft = fittype('1-(x^2)/(a^2)');
% f = fit(rx(1:15)',Rx(1:15)', ft,'StartPoint',[50]);
% l_ux = f.a;
% ft = fittype('1-(x^2)/(a^2)');
% f = fit(ry(1:15)',Ry(1:15)', ft,'StartPoint',[50]);
% l_uy = f.a;
% 
% [Rx,rx,Ry,ry] = chris_autocorrelation(v,v);
% 
% ft = fittype('1-(x^2)/(a^2)');
% f = fit(rx(1:15)',Rx(1:15)', ft,'StartPoint',[50]);
% l_vx = f.a;
% ft = fittype('1-(x^2)/(a^2)');
% f = fit(ry(1:15)',Ry(1:15)', ft,'StartPoint',[50]);
% l_vy = f.a;
% 
% [Rx,rx,Ry,ry] = chris_autocorrelation(C,C);
% 
% ft = fittype('1-(x^2)/(a^2)');
% f = fit(rx(1:15)',Rx(1:15)', ft,'StartPoint',[50]);
% l_cx = f.a;
% ft = fittype('1-(x^2)/(a^2)');
% f = fit(ry(1:15)',Ry(1:15)', ft,'StartPoint',[50]);
% l_cy = f.a;

[dudx,dudy] = chris_gradient(medfilt2(u,[3 3]),1/dpxdx,1/dpxdx);
l_ux2 = (1/dpxdx)*trapz(smooth(sqrt(mean((u.^2)')./mean((dudx.^2)'))));
l_uy2 = (1/dpxdx)*trapz(smooth(sqrt(mean((u.^2)')./mean((dudy.^2)'))));

[dvdx,dvdy] = chris_gradient(medfilt2(v,[3 3]),1/dpxdx,1/dpxdx);
l_vx2 = (1/dpxdx)*trapz(smooth(sqrt(mean((v.^2)')./mean((dvdx.^2)'))));
l_vy2 = (1/dpxdx)*trapz(smooth(sqrt(mean((v.^2)')./mean((dvdy.^2)'))));

[dcdx,dcdy] = chris_gradient(medfilt2(C,[3 3]),1/dpxdx,1/dpxdx);
l_cx2 = (1/dpxdx)*trapz(smooth(sqrt(mean((C.^2)')./mean((dcdx.^2)'))));
l_cy2 = (1/dpxdx)*trapz(smooth(sqrt(mean((C.^2)')./mean((dcdy.^2)'))));

l_cvx2 = (1/dpxdx)*trapz(smooth(sqrt(abs(mean((C.*v)')./mean((dcdx.*dvdx)')))));
l_cvy2 = (1/dpxdx)*trapz(smooth(sqrt(abs(mean((C.*v)')./mean((dcdy.*dvdy)')))));

l_cux2 = (1/dpxdx)*trapz(smooth(sqrt(abs(mean((C.*u)')./mean((dcdx.*dudx)')))));
l_cuy2 = (1/dpxdx)*trapz(smooth(sqrt(abs(mean((C.*u)')./mean((dcdy.*dudy)')))));

l_K2 = (1/dpxdx)*trapz(smooth(sqrt(mean((u.^2)'+(v.^2)')./mean((dudx.^2)'+(dvdy.^2)'))));
l_C2 = (1/dpxdx)*trapz(smooth(sqrt(mean((C.^2)')./mean((dcdx.^2)'+(dcdy.^2)'))));

% l_ux = min([l_ux,l_ux2]);
% l_uy = min([l_uy,l_uy2]);
% 
% l_vx = min([l_vx,l_vx2]);
% l_vy = min([l_vy,l_vy2]);
% 
% l_cx = min([l_cx,l_cx2]);
% l_cy = min([l_cy,l_cy2]);

end

