function [d12,d21] = surfd(img1,img2,searchWin,plot)

% searchWin = 43;
win = (searchWin+1)/2 - 1;

% Binarize Image
BW1 = logical(img1);
BW2 = logical(img2);

% Edge Finding
BW1 = edge(BW1,'sobel');
BW2 = edge(BW2,'sobel');
if (plot == 1)
    figure;
    montage({BW1,BW2})
    title("left: img1 , Right: img2")
end

% d(bw1,bw2)
[x,y] = find(BW1);
d12 = zeros(length(x),1);
for i = 1:length(x)
%     tmp = zeros(size(BW2));
%     tmp(x(i)-win:x(i)+win,y(i)-win:y(i)+win) = BW2(x(i)-win:x(i)+win,y(i)-win:y(i)+win);
    [xTmp,yTmp] = find(BW2);
    dist = [];
    for j = 1:length(xTmp)
        distTmp = sqrt((x(i) - xTmp(j))^2 + (y(i) - yTmp(j))^2);
        dist = [dist, distTmp];
    end
%     min(dist)
    if (~isempty(dist))
        d12(i) = min(dist);
    else
        d12(i) = 9999;
    end
end

% d(bw2,bw1)
[x,y] = find(BW2);
d21 = zeros(length(x),1);
for i = 1:length(x)
%     tmp = zeros(size(BW1));
%     tmp(x(i)-win:x(i)+win,y(i)-win:y(i)+win) = BW1(x(i)-win:x(i)+win,y(i)-win:y(i)+win);
    [xTmp,yTmp] = find(BW1);
    dist = [];
    for j = 1:length(xTmp)
        distTmp = sqrt((x(i) - xTmp(j))^2 + (y(i) - yTmp(j))^2);
        dist = [dist, distTmp];
    end
%     min(dist)
    if (~isempty(dist))
        d21(i) = min(dist);
    else
        d21(i) = 9999;
    end
end
