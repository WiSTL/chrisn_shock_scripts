function [lap,er] = chris_laplacian_2D(u,dx,dy)

[dudx,dudy] = chris_gradient(u,dx,dy);

[d2udx2,d2udxdy1] = chris_gradient(dudx,dx,dy);
[d2udxdy2,d2udy2] = chris_gradient(dudy,dx,dy);

lap = d2udx2+d2udy2;

er = d2udxdy2-d2udxdy1;

end

