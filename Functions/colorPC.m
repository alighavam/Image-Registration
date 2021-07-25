function color = colorPC(mask)

[x,y,z] = ind2sub(size(mask),find(mask));
colorMat = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880 ; 0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840];
color = zeros(length(x),3);
for i = 1:size(color,1)
    tmp = mask(x(i),y(i),z(i));
    color(i,:) = colorMat(tmp,:);
end
color = uint8(color*255);