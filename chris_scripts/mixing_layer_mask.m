function [C1,C2] = mixing_layer_mask(C1,C2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

C1(subplus(C1.*C2)<0.18) = nan;
C2(subplus(C1.*C2)<0.18) = nan;
end

