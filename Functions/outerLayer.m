function [x,y,z,boundary] = outerLayer(mask)

% [~,mask,~] = loadData(fileName);
sz = size(mask);
% mask = flip(mask,3);
[x,y,z] = ind2sub(size(mask),find(mask));
z0 = unique(z);
boundary = zeros(sz);
for i = 1:length(z0)
    inds = find(z == z0(i));
    Xpoints = x(inds);
    Ypoints = y(inds);
    I = zeros(sz(1),sz(2));
    for j = 1:length(Xpoints)
        I(Xpoints(j),Ypoints(j)) = 1;
    end
    BW1 = edge(I,'sobel');
    boundary(:,:,z0(i)) = BW1;
end
[x,y,z] = ind2sub(size(boundary),find(boundary));

