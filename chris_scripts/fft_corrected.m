function [ft1,ft2] = fft_corrected(C,useHann)
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
%     s1(j,:)=(sraw1(j,:)-mean(sraw1(j,:))).*wh;
%     s2(j,:)=(sraw2(j,:)-mean(sraw2(j,:))).*wh;
    s1(j,:)=sraw1(j,:).*wh;
    s2(j,:)=sraw2(j,:).*wh;
end

nfft=2*ncol;
fttemp1=fft(s1)/(wcorrect); % dimension 2 is along row
ft1=fttemp1(:,1:ncol);
% amp1=max(max(abs(ft1)));
fttemp2=fft(s2,nfft,2)/(wcorrect); % dimension 2 is along row
ft2=fttemp2(:,1:ncol)/sqrt(2*ncol*pi^2);
% amp2=max(max(abs(ft2)));
end

