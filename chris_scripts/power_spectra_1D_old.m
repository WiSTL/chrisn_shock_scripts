function [Sc1D,k_x,Sc1Dx] = power_spectra_1D_old(C,dx,antialiase)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

[nx,ny] = size(C);

k_x = (0:nx-1)/dx/nx;

if antialiase

    for i = 1:ny
    
    Sc1Dx_t = (1/nx/dx)*fft(C(:,i));
   
    Sc1Dx_t = Sc1Dx_t(2:floor(end/2));
   
    Sc1D1 = [0 Sc1Dx_t'];
    Sc1D2 = [Sc1Dx_t' 0];
   
    Sc1Dx(:,i) = Sc1D1.*conj(Sc1D2); 
    
    end
else
    
    for i = 1:ny
    
    Sc1Dx = (1/nx/dx)*fft(C(:,i));
   
    Sc1Dx = Sc1Dx(2:floor(end/2));
   
    Sc1Dx(:,i) = Sc1Dx.*conj(Sc1Dx); 
    
    end
    
end

if ny>1
Sc1D = mean(Sc1Dx');
else
   Sc1D =  Sc1Dx;
end
% 
% 
% Sc1Dx = Sc1Dx/ny;
% 
% % Sc1D = Sc1Dx.*conj(Sc1Dx);
% 
% Sc1D = Sc1Dx(2:end/2);
% 
% Sc1D1 = [0 Sc1D'];
% Sc1D2 = [Sc1D' 0];
% 
% Sc1D = Sc1D1.*conj(Sc1D2);
k_x = k_x(2:floor(end/2));

if antialiase
k_x = [k_x 0];
end


end

