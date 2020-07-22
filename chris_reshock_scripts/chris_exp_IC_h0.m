function [h0] = chris_exp_IC_h0(C)

a = mean(C,3);
b = a(500:850,200:900);
c = mean(b');

c = c-min(c);
c = c/(mean(c(1:100)));

[ny,nx] = size(b);

y = [1:ny];

ft = fittype('0.5*erfc((4/sqrt(2*pi))*(x-b)/a)');
f = fit(y(:),c(:), ft,'StartPoint',[200,90]);

h0 = f.a;
y0 = f.b;

d = 0.5*erfc((4/sqrt(2*pi))*(y-y0)/h0);


figure(1)
plot(y,c);
hold on
plot(y,d,'--')
hold off

end

