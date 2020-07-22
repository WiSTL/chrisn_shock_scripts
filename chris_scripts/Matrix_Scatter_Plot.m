function [S,AX,BigAx,H,HAx] = Matrix_Scatter_Plot(M)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[ndp, nf] = size(M);

[S,AX,BigAx,H,HAx]=plotmatrix(M);

n=1;

for i= 1:nf
    
    delete(S(nf*(i-1) + n:i*nf));
    n=n+1;
    
end

for i= 1:nf*nf
    AX(i).Visible = 'off';    
end

for i= 1:nf
    HAx(i).Visible = 'off';    
end


delete(H(1:end));
end

