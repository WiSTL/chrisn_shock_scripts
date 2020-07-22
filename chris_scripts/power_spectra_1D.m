function [espec_ret,k_ret] = power_spectra_1D(C,dx,interlace,useHann)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

sraw1 = C;
sraw2 = sraw1;

[nx,ny]=size(sraw1); % nx number of rows, ny number of columns
nrow=nx;ncol=ny;
% x=dx/2:dx:(nx-1/2)*dx;

i=1:ncol;
wh=hann(ncol)'*(0+useHann)+(1-useHann);   % hanning window, so start and end at same place
wcorrect=sum(wh)/(ncol);
% variable s is for signal, it could be a concentration or velocity
% this loop subtracts the row average and multiplies by hanning window
s1=zeros(nx,ny);
s2=zeros(nx,ny);
for j=1:nrow
    s1(j,:)=(sraw1(j,:)-mean(sraw1(j,:))).*wh;
    s2(j,:)=(sraw2(j,:)-mean(sraw2(j,:))).*wh;
%     s1(j,:)=sraw1(j,:).*wh;
%     s2(j,:)=sraw2(j,:).*wh;
end

nfft=2*ncol;
fttemp1=fft(s1,nfft,2)/(wcorrect); % dimension 2 is along row
ft1=fttemp1(:,1:ncol);
% amp1=max(max(abs(ft1)));
fttemp2=fft(s2,nfft,2)/(wcorrect); % dimension 2 is along row
ft2=fttemp2(:,1:ncol);
% amp2=max(max(abs(ft2)));

k=(2*pi/dx)*(i)/(ncol);
% interlaced energy spectrum, good for experimental data where there is
% noise but unnecessary for simulation data
if nrow > 1
    if interlace == 1
        especi=zeros(nx-1,ny);
        for j=1:nrow-1
            especi(j,:)=ft1(j,:).*conj(ft2(j+1,:)); 
        end;
    else
        especi=zeros(nx,ny);
        for j=1:nrow
            especi(j,:)=ft1(j,:).*conj(ft2(j,:)); 
        end;
    end;
else
        especi(1,:)=ft1(1,:).*conj(ft2(1,:)); 
end
% especmean=mean(especi,1)*ncol/4; % mean interlaced energy spectrum
% espec=especmean;
especi=especi*dx/(2*pi*2*ncol); % mean interlaced energy spectrum
espec_ret=especi;
k_ret=k;

end

