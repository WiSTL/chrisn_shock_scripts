
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Generate Input Images
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 10;
ny = 10;
x = 1:nx;
y = 1:ny;

I=[];
I2=[];

[Y,X] = meshgrid(y,x);

I(:,:,1) = X.^2 + Y.^2;
I(:,:,1) = I(:,:,1)/max(I(:));

I2(:,:,1) = X.^3 - Y.^2;
I2(:,:,1) = I2(:,:,1)/max(I2(:));

I(:,:,2) = squeeze(I(:,:,1))+X.^3;
I(:,:,2) = I(:,:,2)/max(I(:));

I(:,:,3) = squeeze(I(:,:,1))+Y.^4;
I(:,:,3) = I(:,:,3)/max(I(:));

I(:,:,4) = squeeze(I(:,:,1))+(Y.*X).^2;
I(:,:,4) = I(:,:,4)/max(I(:));

GT = 0*X+0.4; %Ground Truth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Genetic Algorithm Inputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = 10; %number of nodes
nl = 5; %number of layers

nb = 20; %number of bots
ng = 100; %number of generations

mp = 0.8;   % mutation probability
pa = 0.2;   % max mutation percent change

score = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Generate Initial Weights
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nb
    [W,B] = chris_GenerateRandomWeights(nn,nl);
    bot{i}.W = W;
    bot{i}.B = B;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   train Neural Net
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[bot,score] = chris_TrainNeuralNet(I,bot,ng,nb,mp,pa,GT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   test Neural Net
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[TopO,TopScore,TopBot] = chris_TestNeuralNet(I2,bot{1},GT);


