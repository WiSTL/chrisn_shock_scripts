function [C,eta,tau] = chris_exp_generate_hs_ensemble()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% read idividual mat files and compose ensemble field 
% in dimensionless coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

eta_interp = -2:0.01:2;
hh_interp = 0:0.1:2;

tau_interp = 0:0.05:5;

dpxdx = 2274.375*2;

n=0;
C_ensemble = [];
for i = 1:fileno
    i
    
     if fileno>1
        file = fileLocations{i};
    else
        file = fileLocations;
     end
   
     load(file);
     
     [C_interp] = chris_interp_hs(C,eta,tau,eta_interp,tau_interp);
     
     for i = 1:nx
         n=n+1;
         C_ensemble(:,n,:) = C_interp(:,i,:);
     end
          
end

C = C_ensemble;
eta = eta_interp;
tau = tau_interp;
end

