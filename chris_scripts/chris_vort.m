function [w,div] = chris_vort(u,v,dx)


[dudy,dudx] = chris_gradient(u,dx,dx);
[dvdy,dvdx] = chris_gradient(v,dx,dx);

w = dvdx-dudy;
div = dudx+dvdy;

end

