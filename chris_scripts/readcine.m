function [ I, H ] = readcine( fname )
%READCINE reads a *.cine and returns a 3D array and a data structure 
% representing the file header. 

%Magic-Numbers----------------
MAGIC = 18755;
%-----------------------------

%Open the cine file
fid = fopen (fname, 'r');
H.Magic = fread(fid, 1, '*uint16');
if (H.Magic ~= MAGIC)
    error(['File "' fname '" is not a Phantom CINE file!']);
end
    
%Read header ----------------------------------------------------------
H.HeadSize = fread (fid, 1, '*uint16');         %Headersize
H.Compression = fread (fid, 1, '*uint16');      %Compression
H.Version = fread (fid, 1, '*uint16');          %Version
H.FirstImageIndex = fread (fid, 1, '*int32');   %FirstMovieImage
H.TotalImages = fread (fid, 1, 'uint32');       %TotalImageCount

firstimageno = fread (fid, 1, 'int32');   	    %FirstImageNo
H.FirstImageNum = firstimageno;

imcount = fread (fid, 1, 'uint32');             %ImageCount
H.ImCount = imcount;

offimheader = fread (fid, 1, 'uint32');         %OffImageHeader
H.OffImageHeader=offimheader;

offsetup = fread (fid, 1, 'uint32');  	  	    %OffSetup
H.OffSetup = offsetup;

offimoffsets = fread (fid, 1, 'uint32');	    %OffImageOffsetscl
H.OffImOffsetsc1 = offimoffsets;

trigfrac = fread (fid, 1, 'uint32');            %Trigger Timing (frac)
H.TriggerFrac = trigfrac;

trigsec = fread (fid, 1, 'uint32');             %Trigger Timing (sec)
H.TriggerSec = trigsec;

fseek (fid, offimheader, 'bof');
H.ImHeadSize=fread (fid, 1, '*uint32');         %Size

width = fread (fid, 1, 'int32');                %Width
H.Width = width;

height = fread (fid, 1, 'int32');               %Height
H.Height = height;

H.Planes = fread (fid, 1, 'uint16');            %Planes

bitdepth = fread (fid, 1, 'uint16');            %8 or 16 bit image
H.BitDepth = bitdepth;

H.Comp = fread (fid, 1, 'uint32');              %Compression
H.SizeImage = fread (fid, 1, 'uint32');         %SizeImage
H.PxPerMX = fread (fid, 1, 'uint32');           %XPelsPerMeter
H.PxPerMY = fread (fid, 1, 'uint32');           %YPelsPerMeter
H.ClrUsed = fread (fid, 1, 'uint32');           %ClrUsed
H.ClrImportant = fread (fid, 1, 'uint32');      %ClrImportant

framerate = fread (fid, 1, 'uint16');           %framerate (fps)
H.FPS = framerate;
%----------------------------------------------------------------------

%Verify read worked
if (imcount < 1)
    error('No images exist in file!');
end 

%Set frame count for array size and get image location
numframes=imcount;
fseek (fid, offimoffsets, 'bof');
imlocs = fread (fid, numframes, 'int64');

%Set bit depth
if (bitdepth == uint16(8))
    bits='uint8';
else
    bits='uint16';
end

%Get location of images
fseek(fid,0,'bof');
for i=1:10000000
    if(fread(fid,1,'*uint16') == 1002)
        break;
    end
end
H.ImageOffset = i;
fseek(fid,2,0);

%Calculate dt for each frame
dt = zeros(numframes,1);
for i=1:numframes
    fracstart = fread (fid, 1, 'uint32');
    secstart = fread (fid, 1, 'uint32');
    dt(i) = (secstart - trigsec) + ((fracstart - trigfrac)/(2^32));
end
H.DeltaT = dt;

%Import actual frames
I=zeros(width,height,numframes);
for i=1:numframes
    fseek (fid, imlocs(i), 'bof');
    jump = fread (fid, 1, 'uint32') - 4;
    fseek (fid, jump, 'cof');
    I(:,:,i) = reshape(fread(fid, width*height, bits), width, height);
end

%Close file
fclose(fid);

end
