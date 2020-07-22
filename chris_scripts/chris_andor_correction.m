function [imgFinal,X] = chris_andor_correction(img)


addpath('img_correction_functions') % Add scripts in functions folder to search path
baseDir = pwd;

%% Settings
frames = 'all';       % Vector of frames to load from image, optionally 'all'
viewFrame = 1;      % Img frame to use for parameter generation
do_bgsub = false;    % Perform background subtraction
Xlimits = [0.05, 0.95]; % Relative acetone concentration lim for img crop
pxcm = 55.5;        % Pixels per cm in imaging plane of input images
pad_px = 1024;      % Number of pixels to pad above/below cropped image

% Set directory containing image files via GUI. 
% Comment out & set directory directly if desired.
imgDir = uigetdir(baseDir, 'Select folder containing images');

% Set to '<image name>.tif to process image other than the one used to
% generate 'params.mat' (located in same imgDir, with same parameters)
imgName = '';

%% Run script
% Attempt to load pre-existing image parameters. 
% Loads 'pxcm', 'imgpath', 'const_molreg', 'origin'
% If do_bgsub = true, also loads 'bgpath', 'bkgd'
%
% 'pxcm': pixels/cm attenuation coef.
% 'imgpath', 'bgpath': absolute paths to image and background files
% 'const_molreg', 'bkgd'::[xmin ymin xmax ymax]: ROIs within image
% 'origin': [x, y] pixel location of laser sheet origin

param_path = setparams(imgDir, do_bgsub,img(:,:,1));
load(param_path);
if ~isempty(imgName)
    imgpath = [imgDir imgName];
end

%% Load image (with optional background subtraction)
% img = tifread(imgpath, 'all');
% bg = tifread(bgpath, 'all');



[ni, nj, n] = size(img); % Number of input images to process

if do_bgsub
    img = bgsub(img,bg);
end

%% Transform image so laser striations are vertical
imgT = lasertransform(img, origin, 'forward');

%% Correct for laser attenuation & non-uniform beam profile
const_y = const_molreg(2):const_molreg(4);
imgCorr = imgT;
A = zeros(n,1);

for i = 1:n
    [imgCorr(:,:,i), A(i)] = beamcorrection(imgT(:,:,i), const_y);
    imgCorr(:,:,i) = imadjust(imgCorr(:,:,i));
    img(:,:,i) = imadjust(img(:,:,i));
end

%% Transform back to Cartesian frame
imgFinal = lasertransform(imgCorr, origin, 'inverse');

%% Find acetone concentration
X{n} = imgFinal(:,:,1); % Initialize output cell array for concentration
for i = 1:n
    X{i} = acetoneX(imgFinal(:,:,i), A(i), pxcm, Xlimits, pad_px);
end


end

