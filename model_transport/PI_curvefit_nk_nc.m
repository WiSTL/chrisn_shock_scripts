k = 0:0.1:100;
lk = log(k);

nc = 11/3;

nk = 5/3;

rel = [100:100:10000];
scl = [0.01:0.1:10];

eLc = 1./sqrt(rel.*scl);%[0.005:0.005:0.06];
eLk = 1./sqrt(rel);%[0.005:0.005:0.05];

[ELC,ELK] = meshgrid(eLc,eLk);


ft = fittype('a*erfc(b*(x-c))');

ac=[];
bc=[];
cc=[];

ak=[];
bk=[];
ck=[];

nn = 0;
for i = eLc
    mm=0;
    nn=nn+1;
    for j = eLk
        mm=mm+1;
        [~,~,pi_c,pi_k] = model_transport(k,nc,nk,i,j);
        pi_c(isnan(pi_c))=0;
        pi_c(isinf(pi_c))=0;
        
        pi_k(isnan(pi_k))=0;
        pi_k(isinf(pi_k))=0;
        
        f = fit(lk(2:end)',pi_c(2:end),ft,'StartPoint',[max(pi_c(2:end)),-1,1]);
        ac(nn,mm) = f.a;
        bc(nn,mm) = f.b;
        cc(nn,mm) = f.c;
        
        f = fit(lk(2:end)',pi_k(2:end),ft,'StartPoint',[max(pi_k(2:end)),-1,1]);
        ak(nn,mm) = f.a;
        bk(nn,mm) = f.b;
        ck(nn,mm) = f.c;
        
    end
end

figure
subplot(3,2,1)
pcolor(scl,rel,ak')
set(gca,'xscale','log')
set(gca,'yscale','log')
axis([10^(-1) 10 ^1 10^2 10^4])
axis square
cmocean('thermal')
shading interp
colorbar

subplot(3,2,2)
pcolor(scl,rel,ac')
set(gca,'xscale','log')
set(gca,'yscale','log')
axis([10^(-1) 10 ^1 10^2 10^4])
axis square
cmocean('thermal')
shading interp
colorbar

subplot(3,2,3)
pcolor(scl,rel,bk')
set(gca,'xscale','log')
set(gca,'yscale','log')
axis([10^(-1) 10 ^1 10^2 10^4])
axis square
cmocean('thermal')
shading interp
colorbar

subplot(3,2,4)
pcolor(scl,rel,bc')
set(gca,'xscale','log')
set(gca,'yscale','log')
axis([10^(-1) 10 ^1 10^2 10^4])
axis square
cmocean('thermal')
shading interp
colorbar

subplot(3,2,5)
pcolor(scl,rel,ck')
set(gca,'xscale','log')
set(gca,'yscale','log')
axis([10^(-1) 10 ^1 10^2 10^4])
axis square
cmocean('thermal')
shading interp
colorbar

subplot(3,2,6)
pcolor(scl,rel,cc')
set(gca,'xscale','log')
set(gca,'yscale','log')
axis([10^(-1) 10 ^1 10^2 10^4])
axis square
cmocean('thermal')
shading interp
colorbar



