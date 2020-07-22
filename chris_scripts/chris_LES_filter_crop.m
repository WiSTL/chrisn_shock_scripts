function [A_filtered] = chris_LES_filter_crop(A,df)

[nz,nx,nt] = size(A);

%,A_fluctuations,A_fluctuations_cropped

df = df/sqrt(8*pi);

for i=1:nt
A_filtered(:,:,i) = imgaussfilt(A(:,:,i),df);
end

% for i=1:nimg
% A_fluctuations(:,:,i) = medfilt2(A(:,:,i)-A_filtered(:,:,i),[df df]);
% end
% 
% for i=1:nimg
% A_fluctuations_cropped(:,:,i) = A_fluctuations(df:end-df,df:end-df,i);
% end


end

