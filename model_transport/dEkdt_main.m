
nk = [1:0.1:4];
nc = [1:0.1:4];

rel = 10:20:500;
scl = 0.1:0.05:10;

k = 0:0.1:100;

Re = 10:50:1000;
Sc = 0.06;
dEk=[];
dEc=[];
dnk=[];
dnc=[];
dlnrel=[];
dlnsc=[];

nn=0;
for i = rel
    mm=0;
    nn=nn+1;
    for j = scl
        mm=mm+1;
        jj=0;
        for l = Re
            jj=jj+1;
            [~,~,~,~,~,~,dEk(nn,mm,jj),dEc(nn,mm,jj),dnk(nn,mm,jj),dnc(nn,mm,jj),dlnrel(nn,mm,jj),dlnsc(nn,mm,jj)] = dEkdt(k,i,j,11/3,5/3,l,Sc);
        end
    end
end
