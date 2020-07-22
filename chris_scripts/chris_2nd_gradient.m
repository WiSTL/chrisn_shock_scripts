function [d2udx2,d2udy2,d2udydx] = chris_2nd_gradient(u,dx,dy)

[dudx,dudy] = chris_gradient(u,dx,dy);

[d2udx2,d2udxdy] = chris_gradient(dudx,dx,dy);

[d2udydx,d2udy2] = chris_gradient(dudy,dx,dy);

d2udydx = 0.5*(d2udydx+d2udxdy);

end

