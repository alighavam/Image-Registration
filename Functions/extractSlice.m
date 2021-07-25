function [dataSlice,labelSlice,inds] = extractSlice(mask,data)

[x,y,z] = ind2sub(size(mask),find(mask));
imgIndx = mode(x);
imgIndy = mode(y);
imgIndz = mode(z);
inds = [imgIndx,imgIndy,imgIndz];
sz = size(mask);
labelSlice{1} = reshape(mask(imgIndx,:,:),sz(2),sz(3))';
labelSlice{2} = reshape(mask(:,imgIndy,:),sz(1),sz(3))';
labelSlice{3} = reshape(mask(:,:,imgIndz),sz(1),sz(2))';
dataSlice{1} = reshape(data(imgIndx,:,:),sz(2),sz(3))';
dataSlice{2} = reshape(data(:,imgIndy,:),sz(1),sz(3))';
dataSlice{3} = reshape(data(:,:,imgIndz),sz(1),sz(2))';