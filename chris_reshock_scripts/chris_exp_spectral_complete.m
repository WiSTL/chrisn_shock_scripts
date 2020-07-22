function [K,Z,Pi_w,Pi_u,Pi_c,Pi_z_w,Pi_z_u,Pi_z_c,P_c,P_u,P_w,E_u,E_w,E_c,D_c,D_u,D_w,D_z_c,D_z_u,D_z_w,X_c,X_u,X_w] = chris_exp_spectral_complete(u,w,c,hl,dpxdx,de,reh,sc)

[nz,nx,nt] = size(c);

% dpxdx = 1; 

k = 2*0.25*dpxdx*[0:nx/2-1,nx/2-1:-1:0]/nx;
z = [1:nz]/dpxdx;

dk = abs(k(2)-k(1));

%[K,Z] = meshgrid(k,z);

[Z,K] = meshgrid(z,k);

cm = squeeze(mean(c,2));
um = squeeze(mean(u,2));
wm = squeeze(mean(w,2));

hw=hann(nx);

for n=1:nt
    for i = 1:nx
        cp(:,i,n) = hw(i)*(c(:,i,n) - smooth(cm(:,n),50));
        up(:,i,n) = hw(i)*(u(:,i,n) - smooth(um(:,n),50));
        wp(:,i,n) = hw(i)*(w(:,i,n) - smooth(wm(:,n),50));
    end
end

for j = 1:nx
    for n = 1:nt
        cmm(:,j,n) = smooth(cm(:,n),50); 
        umm(:,j,n) = smooth(um(:,n),50); 
        wmm(:,j,n) = smooth(wm(:,n),50); 
    end
end

for n = 1:nt
    ch(:,:,n) = (1/nx)*fft(squeeze(cp(:,:,n))')';
    uh(:,:,n) = (1/nx)*fft(squeeze(up(:,:,n))')';
    wh(:,:,n) = (1/nx)*fft(squeeze(wp(:,:,n))')';

    uch(:,:,n) = (1/nx)*fft(squeeze(cp(:,:,n).*up(:,:,n))')';
    wch(:,:,n) = (1/nx)*fft(squeeze(cp(:,:,n).*wp(:,:,n))')';

    u2h(:,:,n) = (1/nx)*fft(squeeze(up(:,:,n).^2)')';
    uwh(:,:,n) = (1/nx)*fft(squeeze(up(:,:,n).*wp(:,:,n))')';
    w2h(:,:,n) = (1/nx)*fft(squeeze(wp(:,:,n).^2)')';
end

for i = 1:nt
    [dw2dz(:,:,i),~] = chris_gradient((w2h(:,:,i)),de,de);
    [duwdz(:,:,i),~] = chris_gradient((uwh(:,:,i)),de,de);
    [dwcdz(:,:,i),~] = chris_gradient((wch(:,:,i)),de,de);
    [dpcdz(:,:,i),~] = chris_gradient(squeeze(cmm(:,:,i)),de,de);
    dpcdz(:,:,i) = (wh(:,:,i).*conj(ch(:,:,i))+ch(:,:,i).*conj(wh(:,:,i))).*dpcdz(:,:,i);
    
    [dpudz(:,:,i),~] = chris_gradient(squeeze(umm(:,:,i)),de,de);
    dpudz(:,:,i) = (wh(:,:,i).*conj(uh(:,:,i))+uh(:,:,i).*conj(wh(:,:,i))).*dpudz(:,:,i);
    
    [dpwdz(:,:,i),~] = chris_gradient(squeeze(wmm(:,:,i)),de,de);
    dpwdz(:,:,i) = (wh(:,:,i).*conj(wh(:,:,i))+wh(:,:,i).*conj(wh(:,:,i))).*dpwdz(:,:,i);
end


for n = 1:nt
    T_w(:,:,n) = -1j*hl*K'.*uwh(:,:,n).*conj(wh(:,:,n));
    T_u(:,:,n) = -1j*hl*K'.*u2h(:,:,n).*conj(uh(:,:,n));
    T_c(:,:,n) = -1j*hl*K'.*uch(:,:,n).*conj(ch(:,:,n));
    
    Pi_z_w(:,:,n) = w2h(:,:,i).*conj(wh(:,:,n))+conj(w2h(:,:,i)).*wh(:,:,n);
    Pi_z_u(:,:,n) = uwh(:,:,i).*conj(uh(:,:,n))+conj(uwh(:,:,i)).*uh(:,:,n);
    Pi_z_c(:,:,n) = wch(:,:,i).*conj(ch(:,:,n))+conj(wch(:,:,i)).*ch(:,:,n);
    
    T_w(:,:,n) = medfilt2(T_w(:,:,n)+conj(T_w(:,:,n)));
    T_u(:,:,n) = medfilt2(T_u(:,:,n)+conj(T_u(:,:,n)));
    T_c(:,:,n) = medfilt2(T_c(:,:,n)+conj(T_c(:,:,n)));
    
    P_c(:,:,n) = medfilt2(dpcdz(:,:,n));
    P_u(:,:,n) = medfilt2(dpudz(:,:,n));
    P_w(:,:,n) = medfilt2(dpwdz(:,:,n));
end

for n = 1:nt
    for j = 1:nz
        Pi_u(j,:,n) = dk*cumtrapz(T_u(j,:,n));
        Pi_w(j,:,n) = dk*cumtrapz(T_w(j,:,n));
        Pi_c(j,:,n) = dk*cumtrapz(T_c(j,:,n));
    end
end

E_u = (uh.*conj(uh));
E_w = (wh.*conj(wh));
E_c = (ch.*conj(ch));

% for i = 1:nt
% [a,K] = pspectrum(squeeze(up(:,:,i))','leakage',0.85);
% E_u(:,:,i) = a';
% E_w(:,:,i) = pspectrum(squeeze(wp(:,:,i))','leakage',0.85)';
% E_c(:,:,i) = pspectrum(squeeze(cp(:,:,i))','leakage',0.85)';
% end
% 
% E_u = imresize3(E_u,[nz,nx,nt]);
% E_w = imresize3(E_w,[nz,nx,nt]);
% E_c = imresize3(E_c,[nz,nx,nt]);
% K = imresize(K,[1,nx]);


for i = 1:nt
    [dEcdz(:,:,i)] = chris_gradient(imgaussfilt(E_c(:,:,i)),de,de);
    [d2Ecdz2(:,:,i)] = chris_gradient(imgaussfilt(dEcdz(:,:,i)),de,de);
    
    [dEudz(:,:,i)] = chris_gradient(imgaussfilt(E_u(:,:,i)),de,de);
    [d2Eudz2(:,:,i)] = chris_gradient(imgaussfilt(dEudz(:,:,i)),de,de);
    
    [dEwdz(:,:,i)] = chris_gradient(imgaussfilt(E_w(:,:,i)),de,de);
    [d2Ewdz2(:,:,i)] = chris_gradient(imgaussfilt(dEwdz(:,:,i)),de,de);
    for j = 1:nx
        D_c(:,j,i) = imgaussfilt(squeeze(-1*(hl*K(j).^2).*E_c(:,j,i)))/(reh*sc);
        D_u(:,j,i) = imgaussfilt(squeeze(-1*(hl*K(j).^2).*E_u(:,j,i))./(reh*(1+6.9.*squeeze(cmm(:,j,i)))));
        D_w(:,j,i) = imgaussfilt(squeeze(-1*(hl*K(j).^2).*E_w(:,j,i))./(reh*(1+6.9.*squeeze(cmm(:,j,i)))));
        
        D_z_c(:,j,i) = imgaussfilt(squeeze(d2Ecdz2(:,j,i)))/(reh*sc);
        D_z_u(:,j,i) = imgaussfilt(squeeze(d2Eudz2(:,j,i))./(reh*(1+6.9.*squeeze(cmm(:,j,i)))));
        D_z_w(:,j,i) = imgaussfilt(squeeze(d2Ewdz2(:,j,i))./(reh*(1+6.9.*squeeze(cmm(:,j,i)))));
    end
end

for i = 1:nt
    [dwpdz(:,:,i),~] = chris_gradient(imgaussfilt(wp(:,:,i)),de,de);
    [dupdz(:,:,i),~] = chris_gradient(imgaussfilt(up(:,:,i)),de,de);
    [dcpdz(:,:,i),~] = chris_gradient(imgaussfilt(cp(:,:,i)),de,de);
    
    dch(:,:,i) = (1/nx)*fft(squeeze(dcpdz(:,:,i))')';
    duh(:,:,i) = (1/nx)*fft(squeeze(dupdz(:,:,i))')';
    dwh(:,:,i) = (1/nx)*fft(squeeze(dwpdz(:,:,i))')';
    
    X_c(:,:,i) = medfilt2(dch(:,:,i).*conj(dch(:,:,i)))/reh;
    X_u(:,:,i) = medfilt2(duh(:,:,i).*conj(duh(:,:,i)))/(reh*sc);
    X_w(:,:,i) = medfilt2(dwh(:,:,i).*conj(dwh(:,:,i)))/(reh*sc);
end

% [Ei,Ki] = meshgrid([-2:0.01:2],k);
% 
% for i = 1:nt
%     [E,K] = meshgrid(squeeze(eta(:,i)),k);
%         
%     a = squeeze(E_c(:,:,i));
%     interpolant = scatteredInterpolant(E(:), K(:), a(:),'natural','none');
%     E_c1(:,:,i) = interpolant(Ei,Ki);
% end





end

