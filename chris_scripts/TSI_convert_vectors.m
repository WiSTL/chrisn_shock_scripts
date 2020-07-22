function [u,v,ug,vg,fileLocations] = TSI_convert_vectors(savefn,roi,dx)
[filename,pathname] = uigetfile('*.vec','multiselect','on');

roicols = roi(1);
roirows = roi(2);

minboxsize = dx;

ncols = round(roicols/minboxsize);
nrows = round(roirows/minboxsize);
filtsize=1;

fileno = size(filename);
fileno = fileno(2);

if fileno>1
    for i = 1:fileno
        fileLocations{i} = [pathname filename{i}];
    end
    
else
    
    fileLocations = [pathname filename];
    
 end

for i = 1:fileno
    i
    
     if fileno>1
        file = fileLocations{i};
    else
        file = fileLocations;
     end
     
    fid=fopen(file,'rt');
    noheadfile = [file,'_noheader','.txt'];
    fid2=fopen(noheadfile,'wt');
    id=0;
    a=fgets(fid);
    while(ischar(a))
        id=id+1;
        if id==1
            a=fgets(fid);
            continue
        else
            fprintf(fid2,a);
        end
        a=fgets(fid);
    end
    
    data = dlmread(noheadfile); 
    
    xlinD=linspace(min(data(:,1)),max(data(:,1)),ncols);
    ylinD=linspace(min(data(:,2)),max(data(:,2)),nrows);
    [X2D,Y2D]=meshgrid(xlinD,ylinD);
    zuxD=griddata(data(:,1),data(:,2),data(:,3),X2D,Y2D,'natural');
    zuyD=griddata(data(:,1),data(:,2),data(:,4),X2D,Y2D,'natural');

    if isempty(data(:,5)) ~= true
    %Uncertainties
    zunsD=griddata(data(:,1),data(:,2),data(:,5),X2D,Y2D,'nearest');
    zuneD=griddata(data(:,1),data(:,2),data(:,6),X2D,Y2D,'nearest');
    su(:,:,i)=medfilt2(zunsD,[filtsize filtsize]);
    eu(:,:,i)=medfilt2(zuneD,[filtsize filtsize]);
    end
    
    v(:,:,i)=medfilt2(zuyD,[filtsize filtsize]);
    u(:,:,i)=medfilt2(zuxD,[filtsize filtsize]);
    
    len=length(data(:,5));
    goodu=NaN(len,1);
    goodz=goodu;
    for j=1:len
        if data(j,5)>0
            goodu(j,1)=data(j,3);
            goodz(j,1)=data(j,4);
        else ploop=1;
        end
    end

    ug=griddata(data(:,1),data(:,2),goodu(:,1),X2D,Y2D,'natural');
    vg=griddata(data(:,1),data(:,2),goodz(:,1),X2D,Y2D,'natural');


end

save(strcat(savefn,'.mat'),'X2D','Y2D','u','v','ug','vg');

if isempty(data(:,5)) ~= true
    save(strcat(savefn,'_uncertainty','.mat'),'X2D','Y2D','eu','su');
end

end



