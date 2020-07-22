function play_img_stack(img,dt)


[~,~,nimg] = size(img);

figure

for i = 1:nimg
    
    imagesc(img(:,:,i));
    axis equal
    axis off
    
    colormap gray
    caxis([0 0.66])
    pause(dt);
    
end
end

