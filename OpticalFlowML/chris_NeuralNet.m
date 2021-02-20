function [O] = chris_NeuralNet(I,W,B)

[ny,nx,input_nodes] = size(I);

[nodes,~,layers] = size(W);

H = zeros(ny,nx,layers,nodes);

for i = 1:nodes
    for j = 1:input_nodes
        H(:,:,1,i) = H(:,:,1,i)+W(i,j,1)*I(:,:,j)+B(i,1);
    end
end

H= tanh(H);

for i = 2:layers
    for j = 1:nodes
        for k = 1:nodes
            H(:,:,i,j) = H(:,:,i,j) + W(j,k,i)*H(:,:,i-1,j)+B(k,i); 
        end
    end
    H = tanh(H);
end
O=0*squeeze(I(:,:,1));
for i = 1:nodes
    O = O + W(i,1,layers)*H(:,:,layers,i)+B(i,layers);
end

O = O/nodes;


end

