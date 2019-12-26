function [Colors] = getColors(NumColor)
Colors = zeros(NumColor ,3);
numPerColor = ceil(NumColor^(1/3));
r = (1:1:numPerColor)/numPerColor;
g = (1:1:numPerColor)/numPerColor;
b = (1:1:numPerColor)/numPerColor;
for index = 1:1:NumColor
    Colors(index,:) = [r(ceil(index/(numPerColor*numPerColor))),g(ceil((mod(index,numPerColor*numPerColor)+1)/numPerColor)),b(mod(mod(index,numPerColor*numPerColor),numPerColor)+1)];
end
