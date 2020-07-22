function [up2,upvp,vp2,vp3,upupvp,upvpdvp,vpvpdup,um,vm,dum,dvm,dupvp,apup,apvp] = reynolds_mean_profile(u,v,eta,rho)

[ny,nx] = size(u);

um = mean(u');
vm = mean(v');

alpha = 1./rho;

am = mean(alpha');
rm = mean(rho');

for j = 1:nx
u_p(:,j) = u(:,j) - um';
end

for j = 1:nx
v_p(:,j) = v(:,j) - vm';
end

for j = 1:nx
a_p(:,j) = alpha(:,j) - am';
end

up2 = mean((u_p.^2)');
vp2 = mean((v_p.^2)');
upvp = mean((u_p.*v_p)');
vp3 = mean((v_p.^3)');
upupvp = mean((u_p.*u_p.*v_p)');
apup = mean((a_p.*u_p)');
apvp = mean((a_p.*v_p)');

[dvp,~] = chris_gradient(v_p,abs(eta(2)-eta(1)),abs(eta(2)-eta(1)));
[dup,~] = chris_gradient(u_p,abs(eta(2)-eta(1)),abs(eta(2)-eta(1)));
[dupvp] = chris_derivative(upvp,abs(eta(2)-eta(1)));

upvpdvp = mean((u_p.*v_p.*dvp)');
vpvpdup = mean((v_p.*v_p.*dup)');

[dvm] = chris_derivative(vm,abs(eta(2)-eta(1)));
[dum] = chris_derivative(um,abs(eta(2)-eta(1)));

end

