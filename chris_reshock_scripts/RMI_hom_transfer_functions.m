[T_hat1,E_hat1,v_hat1,u_hat1,C_hat1,uv_hat1,vv_hat1,kdm1,eta_m1] = chris_exp_1D_homo_fourier(1,[-0.8:0.001:0.5]);
[T_hat2,E_hat2,v_hat2,u_hat2,C_hat2,uv_hat2,vv_hat2,kdm2,eta_m2] = chris_exp_1D_homo_fourier(0,[-0.8:0.001:0.5]);
[T_hat3,E_hat3,v_hat3,u_hat3,C_hat3,uv_hat3,vv_hat3,kdm3,eta_m3] = chris_exp_1D_homo_fourier(0,[-0.8:0.001:0.5]);
[T_hat4,E_hat4,v_hat4,u_hat4,C_hat4,uv_hat4,vv_hat4,kdm4,eta_m4] = chris_exp_1D_homo_fourier(0,[-0.8:0.001:0.5]);

E1 = mean(E_hat1,3);
E2 = mean(E_hat2,3);
E3 = mean(E_hat3,3);
E4 = mean(E_hat4,3);

C1 = mean(C_hat1,3);
C2 = mean(C_hat2,3);
C3 = mean(C_hat3,3);
C4 = mean(C_hat4,3);

G_E12 = E2./E1;
G_E23 = E3./E2;
G_E34 = E4./E3;

G_C12 = (conj(C2).*C2)./(conj(C1).*C1);
G_C23 = (conj(C3).*C3)./(conj(C2).*C2);
G_C34 = (conj(C4).*C4)./(conj(C3).*C3);
