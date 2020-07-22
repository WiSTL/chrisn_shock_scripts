function [en] = chris_exp_cwt(need_dpxdx,eta_range)

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

etam = eta_range;
km = 0.01:0.001:0.3;

for i = 1:fileno
    i
    clear e;
    clear e2;
    clear e3;
    
     if fileno>1
        file = fileLocations{i};
    else
        file = fileLocations;
     end
     
     load(file);
     
    dx_les = 10;
    if need_dpxdx
        if i==15 || i==16 || i==17
            dpxdx = 2274.375*2;
            dx_les = 10;
        else
            dpxdx = 2274.375;
            dx_les = 5;
        end
    end
     
     [ny,nx] = size(C);
     [C,eta,x,delta,y0] = chris_C(C);

    for j = 1:nx
        [e(:,:,j),k] = cwt(C(:,j),'morse','TimeBandwidth',30);
    end
    
    em = mean(e,3);
    
    nk = length(k);
    nke = length(km);
    ne = length(etam);
    
    for j = 1:nk
       e2(j,:) = interp1(eta,em(j,:),etam,'pchip',NaN);        
    end

    e3 = imresize(e2,[nke ne]);
    
    en(i,:,:) = e3;
    
end
end

