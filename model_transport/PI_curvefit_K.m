k = 0:0.1:100;
lk = log(k);

nc = 1;%1-4

nk = [1:0.5:5];

rel = [100:100:10000];

eLk = 1./sqrt(rel);%[0.005:0.005:0.05];

ft = fittype('a*erfc(b*(x-c))');

ak=[];
bk=[];
ck=[];

nn = 0;
for i = nk
    mm=0;
    nn=nn+1;
    for j = eLk
        mm=mm+1;
        [~,~,~,pi_k,~,~,lam_k,~,Ek0] = model_transport(k,nc,i,1,j);
        
        pi_k(isnan(pi_k))=0;
        pi_k(isinf(pi_k))=0;
        
        
        f = fit(lk(2:end)',pi_k(2:end),ft,'StartPoint',[max(pi_k(2:end)),-1,1]);
        ak(nn,mm) = f.a;
        bk(nn,mm) = f.b;
        ck(nn,mm) = f.c;
        lamk(nn,mm) = lam_k;
        E0k(nn,mm) = Ek0;
        
    end
end



