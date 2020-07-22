function [up2n,vp2n,uvn,epsn,Rln,ln] = chris_exp_pdfs(need_dpxdx)

[filename,pathname] = uigetfile('*.mat','multiselect','on');

fileno = size(filename);
fileno=fileno(2)

if fileno>1
    for i = 1:fileno
        fileLocations{i} = [pathname filename{i}];
    end
    
else
    
    fileLocations = [pathname filename];
    
end

un=[];
vn=[];
up2n=[];
vp2n=[];
uvn=[];
cn=[];
cpn=[];
tau12n=[];
t1n=[];
t2n=[];
t3n=[];
xin=[];
etan=[];
etapn=[];
Run=[];
Rvn=[];
P1n=[];
T1n=[];
D1n=[];
e1n=[];
P2n=[];
T2n=[];
D2n=[];
e2n=[];
P3n=[];
T3n=[];
D3n=[];
e3n=[];
tauvn = [];
taumn = [];
Rln = [];
epsn = [];
chin = [];
ln=[];

for i = 1:fileno
    i
    
     if fileno>1
        file = fileLocations{i};
    else
        file = fileLocations;
     end
   
     load(file);
     
     if need_dpxdx
        if i==15 || i==16 || i==17
            dpxdx = 2274.375*2;
            dx_les = 10;
        else
            dpxdx = 2274.375;
            dx_les = 5;
        end
    end
     
     [ny,nx] = size(C);
     [C,eta,x,delta,y0] = chris_C(C);
     [rho,mu] = chris_rho(C);

    [xm,etam] = meshgrid(x,eta);
    
    [delta_dot,v,vm] = chris_v(eta,v,C);

    t0 = h_0/hd_0;
    tau = (1/1000)*t_rs/t0;
    tau_vector = 0*eta+tau;
    tau_matrix = 0*etam+tau;
    
    Cmm = 0.5*erfc(etam);
    Cpm = C-Cmm;
    
    D = 0.00000903221;
    
    area_a = medfilt2(v);
    area_b = medfilt2(v);
    area_a(etam>0)=nan;
    area_b(etam<0)=nan;
    deltadot = nanmean(max(area_a'))-nanmean(min(area_b'));
    
    vm = mean2(v);
    v = v-vm;
    u = u;
    
    Cm = mean(C');
    rm = mean(rho');
    
    [P,T,D_1,e,dd,up2,up3,vp2,vp3,cp2,cp3,apvp] = fluc_equation_terms(u,v,C,rho,deltadot,delta,mu./rho,eta);
    
    C_p = C;
    v_p = v;
    u_p = u;
    r_p = rho;
    vm2 = mean(v');
    um2 = mean(u');
    
    for j = 1:nx
        C_p(:,j) = C(:,j) - Cm';
    end
    
    for j = 1:nx
        r_p(:,j) = rho(:,j) - rm';
    end
    
    for j = 1:nx
        v_p(:,j) = v(:,j) - vm2';
    end
    
    for j = 1:nx
        u_p(:,j) = v(:,j) - um2';
    end
    
    [dupdy,dupdx] = chris_gradient(u_p,1/dpxdx,1/dpxdx);
    [dvpdy,dvpdx] = chris_gradient(v_p,1/dpxdx,1/dpxdx);
    [dcpdy,dcpdx] = chris_gradient(C_p,1/dpxdx,1/dpxdx);
    
    lux = sqrt((u_p.^2)./(dupdx.^2));
    luy = sqrt((u_p.^2)./(dupdy.^2));
    lvx = sqrt((v_p.^2)./(dvpdx.^2));
    lvy = sqrt((v_p.^2)./(dvpdy.^2));
    
    lu = sqrt(lux.^2+luy.^2);
    lv = sqrt(lvx.^2+lvy.^2);
    lambda = medfilt2(sqrt(lu.^2+lv.^2),[10 10]);
    
    Rl = medfilt2(sqrt(v_p.^2+u_p.^2).*lambda./nu_m,[10 10]);   
    
    eps = 2*nu_m*(dupdx.^2+dvpdy.^2);
    eps = medfilt2(eps,[10 10]);
    
    chi = 2*D*(dcpdx.^2+dcpdy.^2);
    chi = medfilt2(chi,[10 10]);
    
    h = delta/dpxdx;
    
    C_e = eps*h/delta_dot^3;
    
    C_x = chi*h/(delta_dot*max(nanmean(C_p.^2,2)));
    
%     figure(1)
%     subplot(2,2,1);
%     imagesc(lambda);
%     axis equal
%     axis off
%     colorbar
%     subplot(2,2,2);
%     imagesc(log(Rl));
%     axis equal
%     axis off
%     colorbar
%     subplot(2,2,3);
%     imagesc(log(eps));
%     axis equal
%     axis off
%     colorbar
%     subplot(2,2,4);
%     imagesc(C);
%     axis equal
%     axis off
%     colorbar
    
    etama = etam(etam>-0.5);
    eta2 = etama(etama<0.5);
    
    Ca = C(etam>-0.5);
    C2 = Ca(etama<0.5);
    
    ua = u(etam>-0.5);
    u2 = ua(etama<0.5);
    
    va = v(etam>-0.5);
    v2 = va(etama<0.5);
    
%     C2 = C(i95:i5,:);
%     u2 = u(i95:i5,:);
%     v2 = v(i95:i5,:);
%     eta2 = etam(i95:i5,:);
%     
%     Ru2 = Ru(i95:i5,:);
%     Rv2 = Rv(i95:i5,:);
    
    un = [un u2(:)'];
    up2n = [up2n up2(:)'];
    vn = [vn v2(:)'];
    vp2n = [vp2n vp2(:)'];
    cn = [cn C2(:)'];
    cpn = [cpn Cpm(:)'];
    tauvn = [tauvn tau_vector(:)'];
    taumn = [taumn tau_matrix(:)'];
    etan = [etan eta2(:)'];
    etapn = [etapn eta(:)'];
    
    Rln = [Rln Rl(:)'];
    epsn = [epsn C_e(:)'];
    chin = [chin C_x(:)'];
    ln = [ln lambda(:)'];
    clear eps;
%     P1n = [P1n P1(:)'];
%     T1n = [T1n T1(:)'];
%     D1n = [D1n D1(:)'];
%     e1n = [e2n e2(:)'];
%     P2n = [P2n P2(:)'];
%     T2n = [T2n T2(:)'];
%     D2n = [D2n D2(:)'];
%     e2n = [e2n e2(:)'];
%     P3n = [P3n P3(:)'];
%     T3n = [T3n T3(:)'];
%     D3n = [D3n D3(:)'];
%     e3n = [e3n e3(:)'];
%     tau12n = [tau12n tau12(:)'];
end

% [P,x,y,cov,r] = chris_jpdf_moments(tauvn,vp2n,1000);
% figure(4)
% pcolor(x,y,medfilt2(P,[3 3]))
% shading flat
% title('P(vp2,up2)')
% xlabel('up2')
% ylabel('vp2')
% colormap(flipud(gray))
% % 

[P,x,y,cov,r] = chris_jpdf_moments(log(epsn),log(chin),2000,[-4 -6],[6 1]);
figure(10)
pcolor(x,y,medfilt2(P,[3 3]))
shading flat
title('P(\eta,v)')
xlabel('\eta')
ylabel('v')
colormap(flipud(gray))

% [P,x,y,cov,r] = chris_jpdf_moments(etan,vn,150,[-1 -50],[1 50]);
% figure(2)
% pcolor(x,y,medfilt2(P,[3 3]))
% shading flat
% title('P(\eta,v)')
% xlabel('\eta')
% ylabel('v')
% colormap(flipud(gray))
% % 
% [P,x,y,cov,r] = chris_jpdf_moments(etan,cn,150,[-1 0],[1 1]);
% figure(3)
% contour(y,x,medfilt2(P,[3 3]),20)
% shading flat
% title('P(\eta,\xi)')
% xlabel('$\xi$')
% ylabel('$\eta$')
% colormap(flipud(gray))
% % % 
% % [P,x,y,cov,r] = chris_jpdf_moments(un,vn,150);
% % figure(4)
% % contour(x,y,medfilt2(P,[3 3]),20)
% % shading flat
% % title('P(u,v)')
% % xlabel('u')
% % ylabel('v')
% % colormap(flipud(gray))
% % 
% % [P,x,y,cov,r] = chris_jpdf_moments(etapn,up2n,1000);
% % figure(5)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(u,v)')
% % xlabel('u')
% % ylabel('v')
% % colormap(flipud(gray))
% % % 
% % [P,x,y,cov,r] = chris_jpdf_moments(etapn,vp2n,1000);
% % figure(6)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(u,v)')
% % xlabel('u')
% % ylabel('v')
% % colormap(flipud(gray))
% % % [P,x,y,cov,r] = chris_jpdf_moments(etapn,P1n,1000);
% % % figure(5)
% % % pcolor(x,y,medfilt2(P,[3 3]))
% % % shading flat
% % % title('P(\eta,P1)')
% % % xlabel('\eta')
% % % ylabel('P1')
% % % colormap(flipud(gray))
% % % 
% % % [P,x,y,cov,r] = chris_jpdf_moments(P2n,P1n,1000);
% % % figure(6)
% % % pcolor(x,y,medfilt2(P,[3 3]))
% % % shading flat
% % % title('P(P2,P1)')
% % % xlabel('P2')
% % % ylabel('P1')
% % % colormap(flipud(gray))
% % 
% % [P,x,y,cov,r] = chris_jpdf_moments(cpn,vn,1000);
% % figure(7)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(C'',v)')
% % xlabel('C''')
% % ylabel('v')
% % colormap(flipud(gray))
% % 
% % [P,x,y,cov,r] = chris_jpdf_moments(cn,un,150);
% % figure(8)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(C'',u)')
% % xlabel('C''')
% % ylabel('u')
% % colormap(flipud(gray))
% % % 
% % [P,x,y,cov,r] = chris_jpdf_moments(cn,vn,150);
% % figure(9)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(C'',v'')')
% % xlabel('C''')
% % ylabel('u')
% % colormap(flipud(gray))
% % 
% % [P,x,y,cov,r] = chris_jpdf_moments(etan,cpn.*vn,1000);
% % figure(10)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(C'',u)')
% % xlabel('\eta')
% % ylabel('c''v')
% % colormap(flipud(gray))
% 
% % [P,x,y,cov,r] = chris_jpdf_moments(cpn,tau12n,1000);
% % figure(14)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(C'',u)')
% % xlabel('C''')
% % ylabel('u')
% % colormap(flipud(gray))
% 
% % [P,x,y,cov,r] = chris_jpdf_moments(log(abs(Run)),log(abs(Rvn)),1000);
% % figure(5)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(Ru,Rv)')
% % xlabel('Ru')
% % ylabel('Rv')
% % colormap(flipud(gray))
% 
% % [P,x,y,cov,r] = chris_jpdf_moments(etan,log(abs(Run)),1000);
% % figure(6)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(Ru,Rv)')
% % xlabel('$\eta$')
% % ylabel('Ru')
% % colormap(flipud(gray))
% 
% % [P,x,y,cov,r] = chris_jpdf_moments(etan,cpn,1000);
% % figure(5)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(\eta,\xi'')')
% % xlabel('\eta')
% % ylabel('\xi''')
% % colormap(flipud(gray))
% 
% % [P,x,y,cov,r] = chris_jpdf_moments(un,cpn,1000);
% % figure(6)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(u,\xi'')')
% % xlabel('u')
% % ylabel('\xi''')
% % colormap(flipud(gray))
% % 
% % [P,x,y,cov,r] = chris_jpdf_moments(vn,cpn,1000);
% % figure(7)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(v,\xi'')')
% % xlabel('v')
% % ylabel('\xi''')
% % colormap(flipud(gray))
% % 
% % 
% % [P,x,y,cov,r] = chris_jpdf_moments(etan,t1n,1000);
% % figure(8)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(\eta,\tau_{xx})')
% % xlabel('\eta')
% % ylabel('\tau_{xx}')
% % colormap(flipud(gray))
% % 
% % [P,x,y,cov,r] = chris_jpdf_moments(etan,xin,1000);
% % figure(9)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(\eta,D)')
% % xlabel('\eta')
% % ylabel('D')
% % colormap(flipud(gray))
% % 
% 
% % size(epsn)
% % size(Rln)
% % 
% % [P,x,y,cov,r] = chris_jpdf_moments(epsn,Rln,100);
% % figure(10)
% % pcolor(x,y,medfilt2(P,[3 3]))
% % shading flat
% % title('P(\eta,uv)')
% % xlabel('\eta')
% % ylabel('uv')
% % colormap(flipud(gray))
% 
% [A,x] = chris_pdf(cn,0.001,nan,nan);
% figure(1)
% plot(x,smooth(A,40))
% xlabel('concentration PDF')
% hold on
% 
% [A,x] = chris_pdf(un,0.001,nan,nan);
% figure(4)
% plot(x,smooth(A,100))
% xlabel('u PDF')
% hold on
% 
% [A,x] = chris_pdf(vn,0.001,nan,nan);
% figure(5)
% plot(x,smooth(A,100))
% xlabel('v PDF')
% hold on
% 
% [A,x] = chris_pdf(epsn,0.0005,nan,nan);
% figure(6)
% plot(x,smooth(A,60))
% xlabel('v PDF')
% hold on
% 
% % [A,x] = chris_pdf(Rln,10000,nan,nan);
% % figure(7)
% % plot(x,smooth(A,20))
% % xlabel('v PDF')
% % hold on
% % 
% % [A,x] = chris_pdf(1000*up2n,0.0001,nan,nan);
% % figure(14)
% % plot(x/1000,smooth(A))
% % xlabel('u'' PDF')
% % hold on
% % 
% % [A,x] = chris_pdf(1000*vp2n,0.0001,nan,nan);
% % figure(15)
% % plot(x/1000,smooth(A))
% % xlabel('v'' PDF')
% % hold on
% % 
% % [A,x] = chris_pdf(1000*uvn,0.0001,nan,nan);
% % figure(16)
% % plot(x/1000,smooth(A))
% % xlabel('u''v'' PDF')
% % hold on


end

