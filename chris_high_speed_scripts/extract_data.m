mean_xt = squeeze(nanmean(c3,2));

[ny,nx,nt] = size(c3);

x = [];
y0 = [];
delta = [];
eta = [];
Cm=[];
cp2=[];
cp3 = [];
cp4 = [];
Cp = [];

for n = 1:nt
    Cm(:,n) = mean(squeeze(c3(:,:,n)),2);
    
    [~,eta(:,n),x(:,n),delta(n),y0(n)] = chris_C(squeeze(c3(:,:,n)));
    
    for i = 1:nx
    Cp(:,i,n) = c3(:,i,n) - Cm(:,n);
    end
    
    cp2(:,n) = mean(squeeze(Cp(:,:,n).^2),2);
    cp3(:,n) = mean(squeeze(Cp(:,:,n).^3),2);
    cp4(:,n) = mean(squeeze(Cp(:,:,n).^4),2);
    
end

t = 1:nt;
e = eta(:,1);
[E,T] = meshgrid(e,t);