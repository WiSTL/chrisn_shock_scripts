function IMAGEStoVID(A,range,filename,dt,cmap)

v = VideoWriter(filename);
v.FrameRate = 1/dt;
open(v);

[nx, ny, nimg] = size(A);

imagesc(flipud(squeeze(A(:,:,1)))) 
axis tight manual 
caxis(range)
cmocean(cmap)
set(gca,'nextplot','replacechildren');

for k = 1:nimg 
    imagesc(flipud(squeeze(A(:,:,k)))) 
    axis equal off 
    caxis(range)
    cmocean(cmap)
    frame = getframe;
    writeVideo(v,frame);
end

close(v);

end

