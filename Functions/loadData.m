function [data,mask,sz] = loadData(fileName)

addpath('./Train')
formatName = '.nii';
dataName = append(fileName,formatName);
labelName = append(fileName,'_label',formatName);

% loadding data
data = im2double(niftiread(dataName));
data = data/max(data(:));
mask = niftiread(labelName);
mask(find(mask)) = mask(find(mask)) - 15;

% [patFlag,pat] = isPat(fileName);
% if (patFlag)
%     data = permute(data,[3,1,2]);
%     mask = permute(mask,[3,1,2]);
%     data = flip(data,2);
%     mask = flip(mask,2);
% end
% if (pat == 20 || pat == 21)
%     data = flip(data,3);
%     mask = flip(mask,3);
% end
sz = size(data);
%     mask = flip(mask,3);