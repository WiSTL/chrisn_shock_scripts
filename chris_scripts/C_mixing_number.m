function [ax,ay,ly,lx] = C_mixing_number(A)
[nx, ny] = size(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ay=zeros(1,nx);
ly=zeros(1,nx);

for i = 1:nx
    m=0;
    n=1;
    by = [];
    by(1)=0;
   for j = 1:ny
       if m==0
           if A(i,j)<0.05
               m=1;
               n=n+1;
               ay(i)=ay(i)+1;
               by(n) = j;
               
           end
       else
           if A(i,j)>0.05
               m=0;
               n=n+1;
               ay(i)=ay(i)+1;
               by(n) = j;
           end
       end
   end
   by1 = [by ny];
   by2 = [0 by];
   lyy = by1-by2;
   
    [~,idx] = max(lyy);
    lyy(idx) = [];
   
   ly(i) = mean(lyy);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax=zeros(1,ny);
lx=zeros(1,ny);

for j = 1:ny
    m=0;
    n=1;
    bx = [];
    bx(1)=0;
   for i = 1:nx
       if m==0
           if A(i,j)<0.05
               m=1;
               n=n+1;
               ax(j)=ax(j)+1;
               bx(n) = i;
               
           end
       else
           if A(i,j)>0.05
               m=0;
               n=n+1;
               ax(j)=ax(j)+1;
               bx(n) = i;
           end
       end
   end
   bx1 = [bx nx];
   bx2 = [0 bx];
   lxx = bx1-bx2;
   
    [~,idx] = max(lxx);
    lxx(idx) = [];
   
   lx(j) = mean(lxx);
end

ax = mean(ax);
ay=mean(ay);

lx = mean(lx)*1/47.39;
ly=mean(ly)*1/47.39;