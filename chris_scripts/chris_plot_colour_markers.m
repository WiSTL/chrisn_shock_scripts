function [] = chris_plot_colour_markers(x,y,z,nf,marker, maxz, minz)

n = length(x);

figure(nf)
hold on

if isnan(maxz)
    maxz = max(z);
end

if isnan(minz)
    minz = min(z);
end

% cm = cmocean('thermal')

for i = 1:n
    c = (z(i) - minz)./(maxz-minz);
    plot(x(i),y(i),'MarkerEdgeColor',[c 0.6*c (1-c)],'Marker',marker)
end

end

