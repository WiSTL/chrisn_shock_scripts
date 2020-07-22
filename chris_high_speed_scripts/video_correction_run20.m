
C_corrected = [];

top(1)=1;
bottom(1) = 10;

top(2)=10;
bottom(2) = 100;

top(3)=50;
bottom(3) = 175;

top(4)=50;
bottom(4) = 250;

top(5)=50;
bottom(5) = 300;

top(6)=50;
bottom(6) = 400;

top(7) = 200;
bottom(7) = 400;

top(8) = 200;
bottom(8) = 500;

top(9) = 400;
bottom(9) = 500;

top(10) = 300;
bottom(10) = 500;

top(11) = 100;
bottom(11) = 400;

top(12) = 100;
bottom(12) = 400;

top(13) = 100;
bottom(13) = 400;

top(13:19) = 50;
bottom(13:19) = 300;

top(20:100) = 10;
bottom(20:100) = 150;

top(101:247) = 100;
bottom(101:247) =400;

[origin, lines] = laserorigin_JGO(squeeze(mean(c1,3)));

for i = 1:244
    i
    C_corrected(:,:,i) = chris_c_correction(squeeze(c1(200:975,:,i)),[top(i) bottom(i) 1 415], origin);
    
end