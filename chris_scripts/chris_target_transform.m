function [trans,itrans,Ro] = chris_target_transform(fixed,moving)

imagesc(fixed);
[x,y] = getpts();

fixedpoints = [x y];

imagesc(moving);
[x,y] = getpts();

movingpoints = [x y];

trans = fitgeotrans(movingpoints,fixedpoints,'nonreflectivesimilarity');
itrans = invert(trans);

Ro = imref2d(size(fixed));
end

