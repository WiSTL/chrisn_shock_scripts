function [snr] = chris_snr(I,bg)

avgI = mean2(I);

stdbg = std2(bg);

snr = avgI/stdbg;


end

