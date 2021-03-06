function out = lasertransform(I, origin, direction)
%% out = lasertransform(I, origin, direction)
% I: 2D or 3D image array to normalize with regards to origin
% origin: (x, y) pixel location of laser sheet origin
% direction: 'forward' or 'inverse'
% out: transformed image of same size as I
intX = origin(1); intY = origin(2);
[nx, ny] = size(I);
out = I;
%out = repmat(I,[1,1,nz]);
movingPoints = [0 0; intX*ny/intY ny; intX ny; intX 0];
fixedPoints =  [0 0; 0            ny; intX ny; intX 0];
tform = fitgeotrans(movingPoints, fixedPoints, 'projective');
R = imref2d([nx,ny]); % Reference image; fixes projection dimensions

if strcmp(direction, 'inverse')
    tform = invert(tform);
end

fillVal = min(min(I(:,:)));
out = imwarp(I, tform, 'Interp', 'cubic', 'FillValues', fillVal);
    
return