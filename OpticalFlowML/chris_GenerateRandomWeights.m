function [W,B] = chris_GenerateRandomWeights(nodes,layers)

for i = 1:nodes
    for j = 1:nodes
        for k = 1:layers
            W(i,j,k) = 2*rand() - 1;
            B(i,k) = 2*rand() - 1;
        end
    end
end

end

