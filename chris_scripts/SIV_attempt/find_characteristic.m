function [X,Y] = find_characteristic(C,x0,y0,dpxdx)

[dcy,dcx] = chris_gradient(C,1/dpxdx,1/dpxdx);

s = 0;
ds = 1;

X = [];
Y = [];

X(1) = x0;
Y(1) = y0;

while dcy(


end

