function [O,score,newbot] = chris_TestNeuralNet(I,bot,GT)

[ny,nx] = size(I);

[O] = chris_NeuralNet(I,bot.W,bot.B);
newbot.O = O;
newbot.W = bot.W;
newbot.B = bot.B;
        
score = sum(abs(O(:) - (GT(:))))/(nx*ny);
newbot.score = score;
    
end

