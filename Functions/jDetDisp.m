function jDetDisp(sx,sy,sz)

det_J = myJacobian(sx,sy,sz);
sprintf("Number of Negative Det = %d",length(find(det_J<0)))
sz = size(det_J);
det_J(find(det_J<0)) = -1;
det_J(find(det_J == 1)) = 0.01;
det_J(find(det_J>0.01)) = 1;

[x,y,z] = ind2sub(size(det_J),find(det_J));
pc = pointCloud([x,y,z]);

color = zeros(length(x),3);
for i = 1:length(x)
    if (det_J(x(i),y(i),z(i)) == 1)
        color(i,:) = [0,1,0];
    elseif (det_J(x(i),y(i),z(i)) == -1)
        color(i,:) = [1,0,0];
    else
        color(i,:) = [0,0,0];
    end
end
detPos = find(color(:,2) == 1);
detNeg = find(color(:,1) == 1);
detZero = find(color(:,1) == 0 & color(:,2) == 0 & color(:,3) == 0);

pcPos = pointCloud([x(detPos),y(detPos),z(detPos)]);
pcNeg = pointCloud([x(detNeg),y(detNeg),z(detNeg)]);
pcTot = pointCloud([x([detPos;detNeg]),y([detPos;detNeg]),z([detPos;detNeg])]);

color = uint8(255*color);
colorPos = color(detPos,:);
colorNeg = color(detNeg,:);
colorTot = color([detPos;detNeg],:);

pcPos.Color = colorPos;
pcNeg.Color = colorNeg;
pcTot.Color = colorTot;

figure;
subplot(1,3,1)
pcshow(pcTot)
title("Positive and Negative Determinan Together")
subplot(1,3,2)
pcshow(pcPos)
title("Positive Determinan")
subplot(1,3,3)
pcshow(pcNeg)
title("Negative Determinan")







