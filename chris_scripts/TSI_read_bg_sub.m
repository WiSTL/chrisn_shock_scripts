function [bg_sub_imgA,bg_sub_imgB] = TSI_read_bg_sub(img_path,bg_path,nimg)

string_endA = '.T000.D000.P000.H000.LA.TIF';
string_endB = '.T000.D000.P000.H000.LB.TIF';

for i = 0:nimg-1
    if i<10
        string_mid =  ['00000' num2str(i)];
    elseif i>=10 && i<100
        string_mid =  ['0000' num2str(i)];
    end
    
   pathA = [img_path string_mid string_endA];
   pathB = [img_path string_mid string_endB];
   
   imgA(:,:,i+1) = tifread(pathA,1);   
   imgB(:,:,i+1) = tifread(pathB,1);
    
end

for i = 0:nimg-1
    if i<10
        string_mid =  ['00000' num2str(i)];
    elseif i>=10 && i<100
        string_mid =  ['0000' num2str(i)];
    end

   pathA = [bg_path string_mid string_endA];
   pathB = [bg_path string_mid string_endB];
   
   bgA(:,:,i+1) = tifread(pathA,1);   
   bgB(:,:,i+1) = tifread(pathB,1);
    
end

bg_sub_imgA = bgsub(imgA,bgA);

bg_sub_imgB = bgsub(imgB,bgB);


end