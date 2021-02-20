function [newbots] = chris_GA(bots,scores,mp,pa,method)

[~,nb] = size(bots);

[nodes,~,layers] = size(bots{1}.W);

[top, topind] = sort(scores, method);

tt = max(1, floor(0.2*nb));

 for j = 1:tt
        newbots{j}.W = bots{topind(1,j)}.W;
        newbots{j}.B = bots{topind(1,j)}.B;
 end

for j = (tt + 1):5*tt
        for k = 1:nodes
            for z = 1:nodes
                for m = 1:layers
                    if rand() < mp
                        newbots{j}.W(k,z,m) = bots{topind(1,j)}.W(k,z,m) + 2*pa*rand() - pa;
                        newbots{j}.B(k,m) = bots{topind(1,j)}.B(k,m) + 2*pa*rand() - pa;
                    end
                end 
            end
        end
end

for j = (5*tt + 1):nb
   [W,B] = chris_GenerateRandomWeights(nodes,layers);
    newbots{j}.W = W;
    newbots{j}.B = B; 
end



