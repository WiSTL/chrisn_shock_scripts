
C_corrected = c1;

top(1)=200;
bottom(1) = 483;

top(2)=200;
bottom(2) = 566;

top(3)=200;
bottom(3) = 662;

top(4)=200;
bottom(4) = 600;

top(5)=655;
bottom(5) = 750;

top(6)=552;
bottom(6) = 738;

top(7) = 442;
bottom(7) = 712;

top(8:247) = 200;
bottom(8:247) = 500;

[origin, lines] = laserorigin_JGO(squeeze(mean(c1,3)));

for i = 1:281
    i
    C_corrected(:,:,i) = chris_c_correction(squeeze(c1(:,:,i)),[top(i) bottom(i) 1 255], origin);
    
end