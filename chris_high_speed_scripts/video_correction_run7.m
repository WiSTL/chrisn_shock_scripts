
C_corrected = [];

top(1)=10;
bottom(1) = 100;

top(2)=10;
bottom(2) = 120;

top(3)=50;
bottom(3) = 170;

top(4)=50;
bottom(4) = 220;

top(5)=50;
bottom(5) = 250;

top(6)=50;
bottom(6) = 250;

top(7) = 10;
bottom(7) = 250;

top(8) = 10;
bottom(8) = 250;

top(9) = 10;
bottom(9) = 230;

top(10) = 10;
bottom(10) = 230;

top(11) = 10;
bottom(11) = 230;

top(12) = 10;
bottom(12) = 230;

top(13:29) = 10;
bottom(13:29) = 160;

top(30:100) = 1;
bottom(30:100) = 50;

top(101:120) = 10;
bottom(101:120) = 80;

top(121:130) = 10;
bottom(121:130) = 100;

top(131:150) = 10;
bottom(131:150) = 200;

[origin, lines] = laserorigin_JGO(squeeze(mean(c1,3)));

for i = 1:150
    i
    C_corrected(:,:,i) = chris_c_correction(squeeze(c1(:,:,i)),[top(i) bottom(i) 1 211], origin);
    
end