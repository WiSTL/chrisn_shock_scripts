
m=1:10;
f = factorial(m);

md=3;

Q=[];

n=[1:1/(2*18):17/3];

eL = [0.05:0.001:0.5];
Rel = eL.^(-2);

ii=0;
for i = n
    ii=ii+1;
    jj=0;
    for j = eL
        jj=jj+1;
        
        F = Fm(m,i,j);
        Fd = diff(F./f);
        
        if Fd(md)<0
            Q(ii,jj)=1;
        else
            Q(ii,jj)=0;
        end
    end
end

