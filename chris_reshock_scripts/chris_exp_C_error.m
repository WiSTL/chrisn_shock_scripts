function [dc] = chris_exp_C_error(yo,xo)

[filename,pathname] = uigetfile('*.mat','multiselect','on');

fileno = size(filename);
fileno=fileno(2)

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
   
     load(file);
     
     [ny,nx] = size(C);
     [C,eta,x,delta,y0] = chris_C(C);
     
    [xm,etam] = meshgrid(x,eta);
    
    dc(i) = mean(nanvar(C(yo,xo)));
    
    figure(2)
    imagesc(C)
    
end

end

