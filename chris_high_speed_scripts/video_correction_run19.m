
C_corrected = [];

top(1)=10;
bottom(1) = 200;

top(2)=10;
bottom(2) = 300;

top(3)=50;
bottom(3) = 400;

top(4)=50;
bottom(4) = 500;

top(5)=50;
bottom(5) = 550;

top(6)=50;
bottom(6) = 600;

top(7) = 200;
bottom(7) = 750;

top(8) = 200;
bottom(8) = 750;

top(9) = 734;
bottom(9) = 800;

top(10) = 640;
bottom(10) = 760;

top(11) = 540;
bottom(11) = 760;

top(12) = 440;
bottom(12) = 750;

top(13) = 345;
bottom(13) = 750;

top(13:29) = 100;
bottom(13:29) = 620;

top(30:100) = 100;
bottom(30:100) = 550;

top(101:247) = 100;
bottom(101:247) =800;

[origin, lines] = laserorigin_JGO(squeeze(mean(c1,3)));

for i = 1:244
    i
    C_corrected(:,:,i) = chris_c_correction(squeeze(c1(:,:,i)),[top(i) bottom(i) 1 415], origin);
    
end