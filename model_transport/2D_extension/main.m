

k = 0:0.5:500;
z = -2:0.1:2;

nu = 5/3;
nw = 5/3;
nc = 11/3;

rel = 10^(0):10^0:5*10^2;
Scl = 1;

eta_L_k = 1./sqrt(rel);%0.05:0.005:0.5;
eta_L_c = 1./sqrt(rel.*Scl);%0.05;

R = 2;
Re = 5000;
Sc = 0.6;

dnu=[];
dnw=[];
dnc=[];
dlnrel=[];
dlnsc=[];

for i = 1:length(rel)
    [~,~,~,~,~,~,~,~,~,dnu(:,i),dnw(:,i),dnc(:,i),dlnrel(:,i),dlnsc(:,i)] = model_transport_terms_2D(k,z,nu,eta_L_k(i),nw,eta_L_k(i),nc,eta_L_c(i),R,Re,Sc);
end