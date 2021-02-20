function [bot,score] = chris_TrainNeuralNet(I,bot,ng,nb,mp,pa,GT)

[ny,nx] = size(I);

for j = 1:ng-1
    for i = 1:nb
        [O] = chris_NeuralNet(I,bot{i}.W,bot{i}.B);
        bot{i}.O = O;
    
        score(i) = sum(abs(O(:) - (GT(:))))/(nx*ny);
        bot{i}.score = score(i);
    end
    
    figure(1)
    plot(score)
    set(gca,'yscale','log')
    pause(0.1)    
    
   [bot] = chris_GA(bot,score,mp,pa,'ascend');    
end

for i = 1:nb
        [O] = chris_NeuralNet(I,bot{i}.W,bot{i}.B);
        bot{i}.O = O;
    
        score(i) = sum(abs(O(:) - (0*O(:)+0.2)))/(nx*ny);
        bot{i}.score = score(i);
end

end

