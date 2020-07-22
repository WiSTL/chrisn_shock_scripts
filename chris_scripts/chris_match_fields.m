function [c,u,v,r] = chris_match_fields(c_corrected,u,v,trans,ro)

[nx ny] = size(u);

c = imresize(chris_transform_concentration(c_corrected, trans, ro),[nx ny]);

u = flipud(u);
v = flipud(v);

imagesc(c); caxis([0 1]);
r = getrect();

c = c(r(2):r(2)+r(4),r(1):r(1)+r(3));
u = u(r(2):r(2)+r(4),r(1):r(1)+r(3));
v = v(r(2):r(2)+r(4),r(1):r(1)+r(3));


end


