
C_corrected = c1;

top(1)=200;
bottom(1) = 373;

top(2)=200;
bottom(2) = 446;

top(3)=200;
bottom(3) = 524;

top(4)=200;
bottom(4) = 597;

top(5)=200;
bottom(5) = 675;

top(6)=710;
bottom(6) = 750;

top(7) = 605;
bottom(7) = 736;

top(8) = 481;
bottom(8) = 716;

top(9) = 481;
bottom(9) = 687;

top(10:29) = 100;
bottom(10:29) = 350;

top(30:100) = 100;
bottom(30:100) = 200;

top(101:247) = 100;
bottom(101:247) = 350;

[origin, lines] = laserorigin_JGO(squeeze(mean(c1,3)));

for i = 1:247
    i
    C_corrected(:,:,i) = chris_c_correction(squeeze(c1(:,:,i)),[top(i) bottom(i) 1 255], origin);
    
end