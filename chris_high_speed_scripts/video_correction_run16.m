
C_corrected = [];

top(1:60)=300;
bottom(1:60) = 700;

top(61:247) = 400;
bottom(61:247) = 900;

[origin, lines] = laserorigin_JGO(squeeze(mean(c1,3)));

for i = 1:244
    i
    C_corrected(:,:,i) = chris_c_correction(squeeze(c1(:,:,i)),[top(i) bottom(i) 1 415], origin);
    
end