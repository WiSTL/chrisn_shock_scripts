function [bg_sub_img,img,bg] = Andor_read_bg_sub(img_path,bg_path)

bg = tifread(bg_path,'all');
img = tifread(img_path,'all');

bg_sub_img = bgsub(img,bg);


end

