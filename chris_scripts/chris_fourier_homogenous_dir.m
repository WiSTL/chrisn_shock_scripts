function [T_hat,E_hat,v_hat,u_hat,C_hat,uv_hat,vv_hat,kd] = chris_fourier_homogenous_dir(C,u,v,delta,eta)

[ny,nx] = size(C);

xx=1:nx;

M = length(xx);
kd = (0:M-1)/M;



for i = 1:ny
    C_hat(i,:)=(1/M)*fft(C(i,:));
    u_hat(i,:)=(1/M)*fft(u(i,:));
    v_hat(i,:)=(1/M)*fft(v(i,:));
    uv_hat(i,:)=(1/M)*fft(u(i,:).*v(i,:));
    uu_hat(i,:)=(1/M)*fft(u(i,:).*u(i,:));
    vv_hat(i,:)=(1/M)*fft(v(i,:).*v(i,:));
end

for i = 1:nx
   dvv_hat(:,i) = chris_derivative(vv_hat(:,i),abs(eta(1)-eta(2)));
   duv_hat(:,i) = chris_derivative(uv_hat(:,i),abs(eta(1)-eta(2)));
end

[etam,kdm] = meshgrid(eta,kd);

kdm = kdm';

T_hat = conj(v_hat).*dvv_hat +conj(u_hat).*duv_hat+1j*kdm.*conj(v_hat).*uv_hat+1j*kdm.*conj(u_hat).*uu_hat;

E_hat = 0.5*(conj(u_hat).*u_hat + conj(v_hat).*v_hat);

end

