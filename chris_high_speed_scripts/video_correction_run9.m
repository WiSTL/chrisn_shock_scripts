
C_corrected = [];

top(1)=10;
bottom(1) = 100;

top(2)=10;
bottom(2) = 150;

top(3)=50;
bottom(3) = 150;

top(4)=50;
bottom(4) = 200;

top(5)=50;
bottom(5) = 300;

top(6)=50;
bottom(6) = 400;

top(7) = 200;
bottom(7) = 400;

top(8) = 200;
bottom(8) = 600;

top(9) = 400;
bottom(9) = 600;

top(10) = 400;
bottom(10) = 750;

top(11) = 630;
bottom(11) = 780;

top(12) = 525;
bottom(12) = 780;

top(13:29) = 100;
bottom(13:29) = 640;

top(30:100) = 100;
bottom(30:100) = 590;

top(101:247) = 100;
bottom(101:247) =800;

[origin, lines] = laserorigin_JGO(squeeze(mean(c1,3)));

for i = 1:244
    i
    C_corrected(:,:,i) = chris_c_correction(squeeze(c1(:,:,i)),[top(i) bottom(i) 1 415], origin);
    
end