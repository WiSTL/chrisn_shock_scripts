function [C] = chris_C_josh_correct(C)

%,eta,x,delta,y0,Cm,x1

    [ny,nx] = size(C);
     
    Cm = squeeze(nanmean(C,2));
    y = [1:ny];
    x = [1:nx];
    
%     ft = fittype('a*exp(b*x)*0.5*erfc((4/sqrt(pi))*(x-d)/c)'); %x-b
%     f = fit(y(:),Cm(:), ft,'StartPoint',[1,-0.002,300,600]);

    al = Cm(10:400);
    bl = y(10:400);

    ft = fittype('a*x+b'); %x-b
    f = fit(bl(:),al(:), ft,'StartPoint',[-0.001,1]);%[1,-0.002,300,600]);
    
    x1 = f.a*y + f.b;
    
%     x1 = f.a.*exp(f.b*y);
%     delta = f.c;
%     y0 = f.d;
    
    for j = 1:nx
        C(:,j) = squeeze(C(:,j))./x1';
    end

%     eta = -1*(y-y0)/delta; %(4/sqrt(2*pi))*

end

