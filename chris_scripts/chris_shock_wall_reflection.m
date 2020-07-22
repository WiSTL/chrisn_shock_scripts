function [Mr,Pr] = chris_shock_wall_reflection(Ms)

Mr = 2.15 - 3.846*exp(-1.207*Ms);

[Pr,~,~,~,~] = normal_shock(Mr,1.664);

end

