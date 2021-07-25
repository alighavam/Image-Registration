% MIAP Project 
% Dr. Fatemizadeh
% Ali Ghavampour - 97102293

% ############## Important ############## 
% 1)
% Change the Name of the Healthy Subject to:
% pat0.nii  && pat0_label.nii

% 2)
%%%% Change MATLAB current folder to the project code folder!!!!
%%%% Then run this part!!!!!
addpath('./Train')
addpath('./Functions')
addpath('./Test')


%% Part 1 - Showing some examples 
clear all; clc;
addpath('./Train')
formatName = '.nii';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = 'pat9'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataName = append(fileName,formatName);
labelName = append(fileName,'_label',formatName);

% loadding data
[data,mask,sz] = loadData(fileName);

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

unique(mask(find(mask)))
color = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880 ; 0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840];

% finding most segmented image
[x,y,z] = ind2sub(size(mask),find(mask));
imgIndx = mode(x);
imgIndy = mode(y);
imgIndz = mode(z);

% Axial Cut
% imgIndz = 160;
img1 = reshape(data(:,:,imgIndz),sz(1),sz(2))'*2;
mask1 = reshape(mask(:,:,imgIndz),sz(1),sz(2))';
RGB = label2rgb(mask1,color,[0.6,0.6,0.6]);
figure;
subplot(1,3,1)
c1 = imfuse(img1,RGB,'blend','scaling','joint');
imshow(c1);
title("Frontal Cut")

% Lateral Cut
img1 = reshape(data(imgIndx,:,:),sz(2),sz(3))'*2;
mask1 = reshape(mask(imgIndx,:,:),sz(2),sz(3))';
RGB = label2rgb(mask1,color,[0.6,0.6,0.6]);
subplot(1,3,2)
c2 = imfuse(img1,RGB,'blend','scaling','joint');
imshow(c2);
title("Lateral Cut")

% Frontal Cut
% imgIndy = 190;
img1 = reshape(data(:,imgIndy,:),sz(1),sz(3))'*2;
mask1 = reshape(mask(:,imgIndy,:),sz(1),sz(3))';
RGB = label2rgb(mask1,color,[0.6,0.6,0.6]);
subplot(1,3,3)
c3 = imfuse(img1,RGB,'blend','scaling','joint');
imshow(c3);
title("Axial Cut")
sgtitle(sprintf("%s",fileName))

verteb = extractVertebra(fileName);
figure;
for i = 1:5
    [x,y,z] = ind2sub(size(verteb{i}),find(verteb{i}));
    shp = alphaShape(x,y,z);
    h = plot(shp);
    h.FaceColor = color(verteb{i}(x(1),y(1),z(1)),:);
    hold on
end
title(sprintf("3D Form of L1 to L5 , %s",fileName))



%% Generating Point Cloud
clear all; clc;
addpath('./Train')
formatName = '.nii';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = 'pat9'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dataName = append(fileName,formatName);
labelName = append(fileName,'_label',formatName);

colorMat = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880 ; 0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840];

% loadding data
data = im2double(niftiread(dataName));
data = data/max(data(:));
mask = niftiread(labelName);
mask(find(mask)) = mask(find(mask)) - 15;

sz = size(data);
mask = flip(mask,3);

[x,y,z] = ind2sub(size(mask),find(mask));
% color = repmat(uint8([0,0,255]),length(x),1);
color = zeros(length(x),3);
for i = 1:size(color,1)
    tmp = mask(x(i),y(i),z(i));
    color(i,:) = colorMat(tmp,:);
end
figure;
color = uint8(color*255);
ptCloud = pointCloud([x,y,z]);
ptCloud.Color = color;
pcshow(ptCloud)
title(sprintf("Point Cloud Plot , %s",fileName))


%% Finding outlayer of data
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = 'pat21'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[~,mask,sz] = loadData(fileName);
mask = flip(mask,3);
[x,y,z] = ind2sub(size(mask),find(mask));
figure
subplot(1,2,1)
pc = pointCloud([x,y,z]);
pcshow(pc)
title("Full points")

[x,y,z,boundary] = outerLayer(mask);
pc = pointCloud([x,y,z]);
subplot(1,2,2)
pcshow(pc)
title("Boundary")
sgtitle(fileName)

%% Extract Vertebras
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = 'pat0'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


verteb = extractVertebra(fileName);
for i = 1:5
    tmp = verteb{i};
    [x,y,z] = ind2sub(size(tmp),find(tmp));
    pc = pointCloud([x,y,z]);
    subplot(3,2,i)
    pcshow(pc)
    title(sprintf("L%d",i))
end


%% L1 to L5 Extraction
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = 'pat0'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


verteb = extractVertebra(fileName);
L1L5 = 0;
for i = 1:5
    L1L5 = L1L5 + verteb{i};
end
[x,y,z] = ind2sub(size(L1L5),find(L1L5));
colorMat = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880 ; 0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840];
color = zeros(length(x),3);
for i = 1:size(color,1)
    tmp = L1L5(x(i),y(i),z(i));
    color(i,:) = colorMat(tmp,:);
end
color = uint8(color*255);
pc = pointCloud([x,y,z]);
pc.Color = color;
pcshow(pc)

%% L1 to L5 Whole Registeration - ICP
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = 'pat9'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


verteb = extractVertebra('pat0'); % can't touch this
atlas = 0;
for i = 1:5
    atlas = atlas + verteb{i};
end
[x,y,z,boundary] = outerLayer(atlas);
pcAtlas = pointCloud([x,y,z]);

verteb = extractVertebra(fileName);
moving = 0;
for i = 1:5
    moving = moving + verteb{i};
end
[x,y,z,boundary] = outerLayer(moving);
pcMoving = pointCloud([x,y,z]);

% downsample
pcAtlas = pcdownsample(pcAtlas,'nonuniformGridSample',10);
pcAtlas.Count
pcMoving = pcdownsample(pcMoving,'nonuniformGridSample',10);
pcMoving.Count

figure;
subplot(1,3,2)
pcshow(pcAtlas)
title("Atlas")

subplot(1,3,1)
pcshow(pcMoving)
title(sprintf("Moving , %s",fileName))

% Registeration 
[tform,movingReg] = pcregistericp(pcMoving,pcAtlas);

subplot(1,3,3)
pcshow(movingReg)
title(sprintf("Registered"))

% Transformed segmentation matrix with segment colors
verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end
[x,y,z] = ind2sub(size(mask),find(mask));
colorMat = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880 ; 0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840];
color = zeros(length(x),3);
for i = 1:size(color,1)
    tmp = mask(x(i),y(i),z(i));
    color(i,:) = colorMat(tmp,:);
end
color = uint8(color*255);
points = [x,y,z];
t = tform.T;
R = t(1:3,1:3);
t = t(4,1:3);
tformPoints = points * R + t;
wholeReg = pointCloud([tformPoints(:,1),tformPoints(:,2),tformPoints(:,3)]);
wholeReg.Color = color;
figure;
subplot(1,3,3)
pcshow(wholeReg)
title("Registeration of whole points")

% Segmentation with segment colors
pcSegment = pointCloud([x,y,z]);
pcSegment.Color = color;
subplot(1,3,1)
pcshow(pcSegment)
title("Original Segmentation")

% Atlas with segment colors
[x,y,z] = ind2sub(size(atlas),find(atlas));
pcAtlas = pointCloud([x,y,z]);
color = zeros(length(x),3);
for i = 1:size(color,1)
    tmp = atlas(x(i),y(i),z(i));
    color(i,:) = colorMat(tmp,:);
end
color = uint8(color*255);
pcAtlas.Color = color;
subplot(1,3,2)
pcshow(pcAtlas);
title("Atlas")

%% Registeration Quality Scores - ICP
clc
% Original Segmenation Mask
verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end
[x,y,z] = ind2sub(size(mask),find(mask));

% Transformed Segmentation Mask
tformedMask = zeros(size(mask));
for i = 1:length(x)
    ind1 = round(tformPoints(i,1));
    ind2 = round(tformPoints(i,2));
    ind3 = round(tformPoints(i,3));
    tformedMask(ind1,ind2,ind3) = mask(x(i),y(i),z(i));
end

% Slicing
[~,sliceAtlas,indsAtlas] = extractSlice(atlas,atlas);
[~,sliceMask,indsMask] = extractSlice_ver2(mask,mask,indsAtlas);
[~,sliceTformedMask,indsTform] = extractSlice_ver2(tformedMask,tformedMask,indsAtlas);


% Scores
% Original Mask and Atlas
searchWin = 43;
dice01 = dice(logical(sliceMask{1}),logical(sliceAtlas{1}));
dice02 = dice(logical(sliceMask{2}),logical(sliceAtlas{2}));
dice03 = dice(logical(sliceMask{3}),logical(sliceAtlas{3}));
hd01 = hd(sliceMask{1},sliceAtlas{1},searchWin);
hd02 = hd(sliceMask{2},sliceAtlas{2},searchWin);
hd03 = hd(sliceMask{3},sliceAtlas{3},searchWin);
asd01 = asd(sliceMask{1},sliceAtlas{1},searchWin);
asd02 = asd(sliceMask{2},sliceAtlas{2},searchWin);
asd03 = asd(sliceMask{3},sliceAtlas{3},searchWin);


% Registered Mask and Atlas
dice11 = dice(logical(sliceTformedMask{1}),logical(sliceAtlas{1}));
dice12 = dice(logical(sliceTformedMask{2}),logical(sliceAtlas{2}));
dice13 = dice(logical(sliceTformedMask{3}),logical(sliceAtlas{3}));
hd11 = hd(sliceTformedMask{1},sliceAtlas{1},searchWin);
hd12 = hd(sliceTformedMask{2},sliceAtlas{2},searchWin);
hd13 = hd(sliceTformedMask{3},sliceAtlas{3},searchWin);
asd11 = asd(sliceTformedMask{1},sliceAtlas{1},searchWin);
asd12 = asd(sliceTformedMask{2},sliceAtlas{2},searchWin);
asd13 = asd(sliceTformedMask{3},sliceAtlas{3},searchWin);

% Dice
fprintf("Dice , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",...
    dice01,dice02,dice03)
fprintf("Dice , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",...
    dice11,dice12,dice13)

% Huasdorff
fprintf("\nHuasdorff , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",...
    hd01,hd02,hd03)
fprintf("Huasdorff , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",...
    hd11,hd12,hd13)

% Huasdorff
fprintf("\nASD , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",...
    asd01,asd02,asd03)
fprintf("ASD , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",...
    asd11,asd12,asd13)

figure;
montage({sliceMask{1},sliceTformedMask{1},sliceAtlas{1},...
    sliceMask{2},sliceTformedMask{2},sliceAtlas{2},...
    sliceMask{3},sliceTformedMask{3},sliceAtlas{3}},'size',[3,3])
title("left: Mask - Middle: ICP Registered - Right: Atlas")
ylabel("Bottom: Frontal , Middle: Axial , Top: Lateral")



%% L1 to L5 Whole Registeration - CPD
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = 'pat9'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


verteb = extractVertebra('pat0'); % Can't touch this
atlas = 0;
for i = 1:5
    atlas = atlas + verteb{i};
end
[x,y,z,boundary] = outerLayer(atlas);

pcAtlas = pointCloud([x,y,z]);

verteb = extractVertebra(fileName);
moving = 0;
for i = 1:5
    moving = moving + verteb{i};
end
[x,y,z,boundary] = outerLayer(moving);
pcMoving = pointCloud([x,y,z]);

% downsample
pcAtlas = pcdownsample(pcAtlas,'nonuniformGridSample',30);
pcAtlas.Count
pcMoving = pcdownsample(pcMoving,'nonuniformGridSample',30);
pcMoving.Count

figure;
subplot(1,3,1)
pcshow(pcAtlas)
title("Atlas")

subplot(1,3,2)
pcshow(pcMoving)
title(sprintf("Moving , %s",fileName))

% Registeration 
[tform,movingReg] = pcregistercpd(pcMoving,pcAtlas,'verbose',true);

subplot(1,3,3)
pcshow(movingReg)
title(sprintf("Registered"))

% Interpolation
clc
[xq,yq,zq] = ind2sub(size(moving),find(moving));
pointsTmp = pcMoving.Location;
x = pointsTmp(:,1);
y = pointsTmp(:,2);
z = pointsTmp(:,3);
regPoints = movingReg.Location;
[x_tformed,y_tformed,z_tformed] = cpdInterpolation(x,y,z,regPoints,xq,yq,zq,'natural');

labelMat = ReconstructInterpMat(x_tformed,y_tformed,z_tformed,xq,yq,zq,moving);
[x,y,z] = ind2sub(size(labelMat),find(labelMat));
color = colorPC(labelMat);
pc = pointCloud([x,y,z]);
pc.Color = color;
figure;
pcshow(pc);
title("Interpolated CPD Transform")


%% Registeration Quality Scores - CPD
clc

% Jacobian
nanInds = find(isnan(x_tformed));
xq(nanInds) = [];
yq(nanInds) = [];
zq(nanInds) = [];
x_tformed(nanInds) = [];
y_tformed(nanInds) = [];
z_tformed(nanInds) = [];
dx = x_tformed - xq;
dy = y_tformed - yq;
dz = z_tformed - zq;
sx = zeros(size(moving));
sy = zeros(size(moving));
sz = zeros(size(moving));
for i = 1:length(dx)
    sx(xq(i),yq(i),zq(i)) = dx(i);
    sy(xq(i),yq(i),zq(i)) = dy(i);
    sz(xq(i),yq(i),zq(i)) = dz(i);
end
jDetDisp(sx,sy,sz);


% Original Segmenation Mask
verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end

% Transformed Segmentation Mask
tformedMask = labelMat;

% Dice Score
[~,sliceAtlas,indsAtlas] = extractSlice(atlas,atlas);
[~,sliceMask,indsMask] = extractSlice_ver2(mask,mask,indsAtlas);
[~,sliceTformedMask,indsTform] = extractSlice_ver2(tformedMask,tformedMask,indsAtlas);

searchWin = 43;
% Original Mask and Atlas
dice01 = dice(logical(sliceMask{1}),logical(sliceAtlas{1}));
dice02 = dice(logical(sliceMask{2}),logical(sliceAtlas{2}));
dice03 = dice(logical(sliceMask{3}),logical(sliceAtlas{3}));
hd01 = hd(sliceMask{1},sliceAtlas{1},searchWin);
hd02 = hd(sliceMask{2},sliceAtlas{2},searchWin);
hd03 = hd(sliceMask{3},sliceAtlas{3},searchWin);
asd01 = asd(sliceMask{1},sliceAtlas{1},searchWin);
asd02 = asd(sliceMask{2},sliceAtlas{2},searchWin);
asd03 = asd(sliceMask{3},sliceAtlas{3},searchWin);

% Registered Mask and Atlas
dice11 = dice(logical(sliceTformedMask{1}),logical(sliceAtlas{1}));
dice12 = dice(logical(sliceTformedMask{2}),logical(sliceAtlas{2}));
dice13 = dice(logical(sliceTformedMask{3}),logical(sliceAtlas{3}));
hd11 = hd(sliceTformedMask{1},sliceAtlas{1},searchWin);
hd12 = hd(sliceTformedMask{2},sliceAtlas{2},searchWin);
hd13 = hd(sliceTformedMask{3},sliceAtlas{3},searchWin);
asd11 = asd(sliceTformedMask{1},sliceAtlas{1},searchWin);
asd12 = asd(sliceTformedMask{2},sliceAtlas{2},searchWin);
asd13 = asd(sliceTformedMask{3},sliceAtlas{3},searchWin);

% Dice
fprintf("Dice , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",dice01,dice02,dice03)
fprintf("Dice , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",dice11,dice12,dice13)

% Huasdorff
fprintf("\nHuasdorff , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",hd01,hd02,hd03)
fprintf("Huasdorff , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",hd11,hd12,hd13)

% Huasdorff
fprintf("\nASD , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",asd01,asd02,asd03)
fprintf("ASD , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",asd11,asd12,asd13)


figure;
montage({sliceMask{1},sliceTformedMask{1},sliceAtlas{1},...
    sliceMask{2},sliceTformedMask{2},sliceAtlas{2},...
    sliceMask{3},sliceTformedMask{3},sliceAtlas{3}},'size',[3,3])
title("left: Mask - Middle: ICP Registered - Right: Atlas")


%% L1 to L5 Whole Registeration - ICP then CPD
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = 'pat9'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


verteb = extractVertebra('pat0'); % Can't touch this 
atlas = 0;
for i = 1:5
    atlas = atlas + verteb{i};
end
[x,y,z,boundary] = outerLayer(atlas);

pcAtlas = pointCloud([x,y,z]);
% pcAtlas.Color = color;

verteb = extractVertebra(fileName);
moving = 0;
for i = 1:5
    moving = moving + verteb{i};
end
[x,y,z,boundary] = outerLayer(moving);
pcMoving = pointCloud([x,y,z]);

% downsample
pcAtlas = pcdownsample(pcAtlas,'nonuniformGridSample',60);
pcAtlas.Count
pcMoving = pcdownsample(pcMoving,'nonuniformGridSample',60);
pcMoving.Count

figure;
subplot(1,3,1)
pcshow(pcAtlas)
title("Atlas")

subplot(1,3,2)
pcshow(pcMoving)
title(sprintf("Moving , %s",fileName))

% Registeration 
[tformICP,movingRegICP] = pcregistericp(pcMoving,pcAtlas);
[tformCPD,movingReg] = pcregistercpd(movingRegICP,pcAtlas,'verbose',true);

subplot(1,3,3)
pcshow(movingReg)
title(sprintf("Registered"))

% ICP Transformation for all points
verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end
[x,y,z] = ind2sub(size(mask),find(mask));
points = [x,y,z];
t = tformICP.T;
R = t(1:3,1:3);
t = t(4,1:3);
tformPoints = points * R + t;

% CPD Interpolation
xq = tformPoints(:,1);
yq = tformPoints(:,2);
zq = tformPoints(:,3);
pointsTmp = movingRegICP.Location;
x = pointsTmp(:,1);
y = pointsTmp(:,2);
z = pointsTmp(:,3);
regPoints = movingReg.Location;
[x_tformed,y_tformed,z_tformed] = cpdInterpolation(x,y,z,regPoints,xq,yq,zq,'natural');

% Visualization on interp CPD
pc = pointCloud([x_tformed,y_tformed,z_tformed]);
figure;
pcshow(pc);
title("Interpolated CPD Transform")


%% Registeration Quality Scores - ICP then CPD
clc

% Jacobian
xq = points(:,1);
yq = points(:,2);
zq = points(:,3);
nanInds = find(isnan(x_tformed));
xq(nanInds) = [];
yq(nanInds) = [];
zq(nanInds) = [];
x_tformed(nanInds) = [];
y_tformed(nanInds) = [];
z_tformed(nanInds) = [];
dx = x_tformed - xq;
dy = y_tformed - yq;
dz = z_tformed - zq;
sx = zeros(size(mask));
sy = zeros(size(mask));
sz = zeros(size(mask));
for i = 1:length(dx)
    sx(xq(i),yq(i),zq(i)) = dx(i);
    sy(xq(i),yq(i),zq(i)) = dy(i);
    sz(xq(i),yq(i),zq(i)) = dz(i);
end
jDetDisp(sx,sy,sz);


tformPoints = pc.Location;
nanInds = find(isnan(tformPoints(:,1)));

% Original Segmenation Mask
verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end
[x,y,z] = ind2sub(size(mask),find(mask));
x(nanInds) = [];
y(nanInds) = [];
z(nanInds) = [];

% Transformed Segmentation Mask
tformedMask = zeros(size(mask));
tformPoints(nanInds,:) = [];
for i = 1:length(x)
    ind1 = round(tformPoints(i,1));
    ind2 = round(tformPoints(i,2));
    ind3 = round(tformPoints(i,3));
    tformedMask(ind1,ind2,ind3) = mask(x(i),y(i),z(i));
end

% Dice Score
[~,sliceAtlas,indsAtlas] = extractSlice(atlas,atlas);
[~,sliceMask,indsMask] = extractSlice_ver2(mask,mask,indsAtlas);
[~,sliceTformedMask,indsTform] = extractSlice_ver2(tformedMask,tformedMask,indsAtlas);

searchWin = 43;
% Original Mask and Atlas
dice01 = dice(logical(sliceMask{1}),logical(sliceAtlas{1}));
dice02 = dice(logical(sliceMask{2}),logical(sliceAtlas{2}));
dice03 = dice(logical(sliceMask{3}),logical(sliceAtlas{3}));
hd01 = hd(sliceMask{1},sliceAtlas{1},searchWin);
hd02 = hd(sliceMask{2},sliceAtlas{2},searchWin);
hd03 = hd(sliceMask{3},sliceAtlas{3},searchWin);
asd01 = asd(sliceMask{1},sliceAtlas{1},searchWin);
asd02 = asd(sliceMask{2},sliceAtlas{2},searchWin);
asd03 = asd(sliceMask{3},sliceAtlas{3},searchWin);

% Registered Mask and Atlas
dice11 = dice(logical(sliceTformedMask{1}),logical(sliceAtlas{1}));
dice12 = dice(logical(sliceTformedMask{2}),logical(sliceAtlas{2}));
dice13 = dice(logical(sliceTformedMask{3}),logical(sliceAtlas{3}));
hd11 = hd(sliceTformedMask{1},sliceAtlas{1},searchWin);
hd12 = hd(sliceTformedMask{2},sliceAtlas{2},searchWin);
hd13 = hd(sliceTformedMask{3},sliceAtlas{3},searchWin);
asd11 = asd(sliceTformedMask{1},sliceAtlas{1},searchWin);
asd12 = asd(sliceTformedMask{2},sliceAtlas{2},searchWin);
asd13 = asd(sliceTformedMask{3},sliceAtlas{3},searchWin);

% Dice
fprintf("Dice , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",dice01,dice02,dice03)
fprintf("Dice , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",dice11,dice12,dice13)

% Huasdorff
fprintf("\nHuasdorff , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",hd01,hd02,hd03)
fprintf("Huasdorff , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",hd11,hd12,hd13)

% Huasdorff
fprintf("\nASD , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",asd01,asd02,asd03)
fprintf("ASD , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",asd11,asd12,asd13)


figure;
montage({sliceMask{1},sliceTformedMask{1},sliceAtlas{1},...
    sliceMask{2},sliceTformedMask{2},sliceAtlas{2},...
    sliceMask{3},sliceTformedMask{3},sliceAtlas{3}},'size',[3,3])
title("left: Mask - Middle: ICP Registered - Right: Atlas")



%% Single Vertebra Registeration - ICP
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = 'pat9'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


atlas = extractVertebra('pat0'); % can't touch this
pat = extractVertebra(fileName);

% Visualization
for i = 1:5
    [x,y,z] = ind2sub(size(atlas{i}),find(atlas{i}));
    pcAtlas{i} = pointCloud([x,y,z]);
    subplot(2,5,i);
    pcshow(pcAtlas{i});
    title(sprintf("Atlas - L%d",i))
    
    [x,y,z] = ind2sub(size(pat{i}),find(pat{i}));
    pcPat{i} = pointCloud([x,y,z]);
    subplot(2,5,i+5)
    pcshow(pcPat{i})
    title(sprintf("%s - L%d",fileName,i))
end

% Extracting Outlayer
for i = 1:5
    [~,~,~,boundary] = outerLayer(atlas{i});
    atlas{i} = boundary;
    [~,~,~,boundary] = outerLayer(pat{i});
    pat{i} = boundary;
end

% ICP Registeration
for i = 1:5
    disp(i)
    [x,y,z] = ind2sub(size(atlas{i}),find(atlas{i}));
    pcTmpFixed = pointCloud([x,y,z]);
%     pcTmpFixed = pcdownsample(pcTmpFixed,'nonuniformGridSample',6);
    pcTmpFixed.Count
    pcTmpFixed.Normal = pcnormals(pcTmpFixed);
    [x,y,z] = ind2sub(size(pat{i}),find(pat{i}));
    pcTmpMoving = pointCloud([x,y,z]);
%     pcTmpMoving = pcdownsample(pcTmpMoving,'nonuniformGridSample',50);
    pcTmpMoving.Count
    pcTmpMoving.Normal = pcnormals(pcTmpMoving);
    [tform{i},pcReg{i}] = pcregistericp(pcTmpMoving,pcTmpFixed);
end

figure;
for i = 1:5
    subplot(1,5,i)
    pcshow(pcReg{i});
    title(sprintf("ICP Registered L%d",i))
end

% Merge single vertebras
pcRegSave = pcReg;
location = [];
for i = 1:5
    pcTmp = pcReg{i};
    location = [location ; pcTmp.Location];
end
pcMerged = pointCloud(location);
figure;
subplot(1,3,3)
pcshow(pcMerged)
title("Merged Point Cloud - ICP")

% Atlas Point Cloud
subplot(1,3,2)
verteb = extractVertebra('pat0');
atlas = 0;
for i = 1:5
    atlas = atlas + verteb{i};
end
[x,y,z] = ind2sub(size(atlas),find(atlas));
pcAtlas = pointCloud([x,y,z]);
pcshow(pcAtlas)
title("Atlas")

% Segment Point Cloud
subplot(1,3,1)
verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end
[x,y,z] = ind2sub(size(mask),find(mask));
pcMask = pointCloud([x,y,z]);
pcshow(pcMask)
title(fileName)

% Interpolation
[xq,yq,zq] = ind2sub(size(mask),find(mask));
x = [];
y = [];
z = [];
for i = 1:5
    [xTmp,yTmp,zTmp] = ind2sub(size(pat{i}),find(pat{i}));
    x = [x ; xTmp];
    y = [y ; yTmp];
    z = [z ; zTmp];
end

regPoints = pcMerged.Location;
[x_tformed,y_tformed,z_tformed] = cpdInterpolation(x,y,z,regPoints,xq,yq,zq,'natural');


% Visualization
labelMat = ReconstructInterpMat(x_tformed,y_tformed,z_tformed,xq,yq,zq,mask);
[x,y,z] = ind2sub(size(labelMat),find(labelMat));
color = colorPC(labelMat);
pcReg = pointCloud([x_tformed,y_tformed,z_tformed]);

pcTmp = pointCloud([x,y,z]);
pcTmp.Color = color;
figure;
pcshow(pcTmp);
title("Interpolated ICP Transform")


%% Registeration Quality Scores - Single Verterba ICP
clc

% Jacobian
nanInds = find(isnan(x_tformed));
xq(nanInds) = [];
yq(nanInds) = [];
zq(nanInds) = [];
x_tformed(nanInds) = [];
y_tformed(nanInds) = [];
z_tformed(nanInds) = [];
dx = x_tformed - xq;
dy = y_tformed - yq;
dz = z_tformed - zq;
sx = zeros(size(mask));
sy = zeros(size(mask));
sz = zeros(size(mask));
for i = 1:length(dx)
    sx(xq(i),yq(i),zq(i)) = dx(i);
    sy(xq(i),yq(i),zq(i)) = dy(i);
    sz(xq(i),yq(i),zq(i)) = dz(i);
end
jDetDisp(sx,sy,sz);


tformPoints = pcReg.Location;
nanInds = find(isnan(tformPoints(:,1)));

% Original Segmenation Mask
verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end
[x,y,z] = ind2sub(size(mask),find(mask));
x(nanInds) = [];
y(nanInds) = [];
z(nanInds) = [];

% Transformed Segmentation Mask
tformedMask = zeros(size(mask));
tformPoints(nanInds,:) = [];
for i = 1:length(x)
    ind1 = round(tformPoints(i,1));
    ind2 = round(tformPoints(i,2));
    ind3 = round(tformPoints(i,3));
    tformedMask(ind1,ind2,ind3) = mask(x(i),y(i),z(i));
end

% Dice Score
[~,sliceAtlas,indsAtlas] = extractSlice(atlas,atlas);
[~,sliceMask,indsMask] = extractSlice_ver2(mask,mask,indsAtlas);
[~,sliceTformedMask,indsTform] = extractSlice_ver2(tformedMask,tformedMask,indsAtlas);

searchWin = 43;
% Original Mask and Atlas
dice01 = dice(logical(sliceMask{1}),logical(sliceAtlas{1}));
dice02 = dice(logical(sliceMask{2}),logical(sliceAtlas{2}));
dice03 = dice(logical(sliceMask{3}),logical(sliceAtlas{3}));
hd01 = hd(sliceMask{1},sliceAtlas{1},searchWin);
hd02 = hd(sliceMask{2},sliceAtlas{2},searchWin);
hd03 = hd(sliceMask{3},sliceAtlas{3},searchWin);
asd01 = asd(sliceMask{1},sliceAtlas{1},searchWin);
asd02 = asd(sliceMask{2},sliceAtlas{2},searchWin);
asd03 = asd(sliceMask{3},sliceAtlas{3},searchWin);

% Registered Mask and Atlas
dice11 = dice(logical(sliceTformedMask{1}),logical(sliceAtlas{1}));
dice12 = dice(logical(sliceTformedMask{2}),logical(sliceAtlas{2}));
dice13 = dice(logical(sliceTformedMask{3}),logical(sliceAtlas{3}));
hd11 = hd(sliceTformedMask{1},sliceAtlas{1},searchWin);
hd12 = hd(sliceTformedMask{2},sliceAtlas{2},searchWin);
hd13 = hd(sliceTformedMask{3},sliceAtlas{3},searchWin);
asd11 = asd(sliceTformedMask{1},sliceAtlas{1},searchWin);
asd12 = asd(sliceTformedMask{2},sliceAtlas{2},searchWin);
asd13 = asd(sliceTformedMask{3},sliceAtlas{3},searchWin);

% Dice
fprintf("Dice , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",dice01,dice02,dice03)
fprintf("Dice , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",dice11,dice12,dice13)

% Huasdorff
fprintf("\nHuasdorff , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",hd01,hd02,hd03)
fprintf("Huasdorff , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",hd11,hd12,hd13)

% Huasdorff
fprintf("\nASD , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",asd01,asd02,asd03)
fprintf("ASD , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",asd11,asd12,asd13)


% intersection of single Vertebras
for i = 1:4
    reg1 = pcRegSave{i};
    reg1 = reg1.Location;
    reg2 = pcRegSave{i+1};
    reg2 = reg2.Location;
    shp1 = alphaShape(reg1(:,1),reg1(:,2),reg1(:,3));
    shp2 = alphaShape(reg2(:,1),reg2(:,2),reg2(:,3));
    id1=inShape(shp2,reg1(:,1),reg1(:,2),reg1(:,3));
    id2=inShape(shp1,reg2(:,1),reg2(:,2),reg2(:,3));
    shp3=alphaShape([reg1(id1,1); reg2(id2,1)], [reg1(id1,2); reg2(id2,2)], ...
        [reg1(id1,3); reg2(id2,3)]);
    shpHolder{i} = shp3;
    V(i) = volume(shp3);
end
fprintf("\nVolumes:\n L1 & L2 = %.4f , L2 & L3 = %.4f , L3 & L4 = %.4f , L4 & L5 = %.4f\n",...
    V(1),V(2),V(3),V(4))

figure;
for i = 1:5
    reg = pcRegSave{i};
    reg = reg.Location;
    shp = alphaShape(reg(:,1),reg(:,2),reg(:,3));
    h = plot(shp);
    h.FaceAlpha = 0.05;
    hold on
end
for i = 1:4
    h = plot(shpHolder{i});
    h.FaceColor = 'r';
    hold on
end
title("Intersections After Single Vertebra ICP")

figure;
montage({sliceMask{1},sliceTformedMask{1},sliceAtlas{1},...
    sliceMask{2},sliceTformedMask{2},sliceAtlas{2},...
    sliceMask{3},sliceTformedMask{3},sliceAtlas{3}},'size',[3,3])
title("left: Mask - Middle: ICP Registered - Right: Atlas")



%% Single Vertebra Registeration - CPD
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = 'pat9'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


atlas = extractVertebra('pat0'); % can't touch this
pat = extractVertebra(fileName);

% Visualization
for i = 1:5
    [x,y,z] = ind2sub(size(atlas{i}),find(atlas{i}));
    pcAtlas{i} = pointCloud([x,y,z]);
    subplot(2,5,i);
    pcshow(pcAtlas{i});
    title(sprintf("Atlas - L%d",i))
    
    [x,y,z] = ind2sub(size(pat{i}),find(pat{i}));
    pcPat{i} = pointCloud([x,y,z]);
    subplot(2,5,i+5)
    pcshow(pcPat{i})
    title(sprintf("%s - L%d",fileName,i))
end

% Extracting Outlayer
for i = 1:5
    [~,~,~,boundary] = outerLayer(atlas{i});
    atlas{i} = boundary;
    [~,~,~,boundary] = outerLayer(pat{i});
    pat{i} = boundary;
end

% CPD Registeration
for i = 1:5
    disp(i)
    [x,y,z] = ind2sub(size(atlas{i}),find(atlas{i}));
    pcTmpFixed = pointCloud([x,y,z]);
    pcTmpFixed = pcdownsample(pcTmpFixed,'nonuniformGridSample',12);
    pcTmpFixed.Count
    [x,y,z] = ind2sub(size(pat{i}),find(pat{i}));
    pcTmpMoving = pointCloud([x,y,z]);
    pcTmpMoving = pcdownsample(pcTmpMoving,'nonuniformGridSample',12);
    pcMovingBasePoints{i} = pcTmpMoving.Location;
    pcTmpMoving.Count
    [tform{i},pcReg{i}] = pcregistercpd(pcTmpMoving,pcTmpFixed,'verbose',true);
end

figure;
for i = 1:5
    subplot(1,5,i)
    pcshow(pcReg{i});
    title(sprintf("CPD Registered L%d",i))
end

% Merge single vertebras
location = [];
for i = 1:5
    pcTmp = pcReg{i};
    location = [location ; pcTmp.Location];
end
pcMerged = pointCloud(location);
figure;
subplot(1,3,3)
pcshow(pcMerged)
title("Merged Registered Vertebras - CPD")

% Atlas Point Cloud
subplot(1,3,2)
verteb = extractVertebra('pat0');
atlas = 0;
for i = 1:5
    atlas = atlas + verteb{i};
end
[x,y,z] = ind2sub(size(atlas),find(atlas));
pcAtlas = pointCloud([x,y,z]);
pcshow(pcAtlas)
title("Atlas")

% Segment Point Cloud
subplot(1,3,1)
verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end
[x,y,z] = ind2sub(size(mask),find(mask));
pcMask = pointCloud([x,y,z]);
pcshow(pcMask)
title(fileName)


% Interpolation of Single Vertebra - CPD
clc;
regPoints = [];
for i = 1:5
    tmp = pcReg{i};
    tmp = tmp.Location;
    regPoints = [regPoints ; tmp];
end

x = [];
y = [];
z = [];
for i = 1:5
    x = [x ; pcMovingBasePoints{i}(:,1)];
    y = [y ; pcMovingBasePoints{i}(:,2)];
    z = [z ; pcMovingBasePoints{i}(:,3)];
end

% Grid Points
verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end
[xq,yq,zq] = ind2sub(size(mask),find(mask));

% interpolation for X dimension --> x_tformed = f1(x,y,z), we will find f1
tic
x_tformed = griddata(x,y,z,regPoints(:,1),xq,yq,zq,'natural');
toc

% interpolation for Y dimension --> y_tformed = f2(x,y,z), we will find f2
tic
y_tformed = griddata(x,y,z,regPoints(:,2),xq,yq,zq,'natural');
toc

% interpolation for Z dimension --> z_tformed = f3(x,y,z), we will find f3
tic
z_tformed = griddata(x,y,z,regPoints(:,3),xq,yq,zq,'natural');
toc

pcAllPoints = pointCloud([x_tformed,y_tformed,z_tformed]);
figure;
pcshow(pcAllPoints);

% Reconstruct Label Matrix From Interpolated Points
labelMat = ReconstructInterpMat(x_tformed,y_tformed,z_tformed,xq,yq,zq,mask);


% Visulization:
% Atlas
verteb = extractVertebra('pat0');
atlas = 0;
for i = 1:5
    atlas = atlas + verteb{i};
end
[x,y,z] = ind2sub(size(atlas),find(atlas));
pcAtlas = pointCloud([x,y,z]);
color = colorPC(atlas);
pcAtlas.Color = color;
figure;
subplot(1,3,2)
pcshow(pcAtlas)
title("Atlas")

% labelMat
[x,y,z] = ind2sub(size(labelMat),find(labelMat));
pcLabelMat = pointCloud([x,y,z]);
color = colorPC(labelMat);
pcLabelMat.Color = color;
subplot(1,3,3)
pcshow(pcLabelMat)
title("CPD Registered With Full Poitns")

% Original Mask
verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end
[x,y,z] = ind2sub(size(mask),find(mask));
pc = pointCloud([x,y,z]);
color = colorPC(mask);
pc.Color = color;
subplot(1,3,1)
pcshow(pc)
title("Oiginal Mask")


%% Registeration Quality Scores - Single Vertebra CPD
clc

tformPoints = [x_tformed,y_tformed,z_tformed];
nanInds = find(isnan(tformPoints(:,1)));

% Jacobian
x_tformedJ = x_tformed;
y_tformedJ = y_tformed;
z_tformedJ = z_tformed;
nanIndsJ = find(isnan(x_tformed));
xq = x;
yq = y;
zq = z;
xq(nanIndsJ) = [];
yq(nanIndsJ) = [];
zq(nanIndsJ) = [];
x_tformedJ(nanIndsJ) = [];
y_tformedJ(nanIndsJ) = [];
z_tformedJ(nanIndsJ) = [];
dx = x_tformedJ - xq;
dy = y_tformedJ - yq;
dz = z_tformedJ - zq;
sx = zeros(size(mask));
sy = zeros(size(mask));
sz = zeros(size(mask));
for i = 1:length(dx)
    sx(xq(i),yq(i),zq(i)) = dx(i);
    sy(xq(i),yq(i),zq(i)) = dy(i);
    sz(xq(i),yq(i),zq(i)) = dz(i);
end
jDetDisp(sx,sy,sz);


% Original Segmenation Mask
verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end
[x,y,z] = ind2sub(size(mask),find(mask));
x(nanInds) = [];
y(nanInds) = [];
z(nanInds) = [];

% Transformed Segmentation Mask
tformedMask = zeros(size(mask));
tformPoints(nanInds,:) = [];
for i = 1:length(x)
    ind1 = round(tformPoints(i,1));
    ind2 = round(tformPoints(i,2));
    ind3 = round(tformPoints(i,3));
    tformedMask(ind1,ind2,ind3) = mask(x(i),y(i),z(i));
end

% Dice Score
[~,sliceAtlas,indsAtlas] = extractSlice(atlas,atlas);
[~,sliceMask,indsMask] = extractSlice_ver2(mask,mask,indsAtlas);
[~,sliceTformedMask,indsTform] = extractSlice_ver2(tformedMask,tformedMask,indsAtlas);

searchWin = 43;
% Original Mask and Atlas
dice01 = dice(logical(sliceMask{1}),logical(sliceAtlas{1}));
dice02 = dice(logical(sliceMask{2}),logical(sliceAtlas{2}));
dice03 = dice(logical(sliceMask{3}),logical(sliceAtlas{3}));
hd01 = hd(sliceMask{1},sliceAtlas{1},searchWin);
hd02 = hd(sliceMask{2},sliceAtlas{2},searchWin);
hd03 = hd(sliceMask{3},sliceAtlas{3},searchWin);
asd01 = asd(sliceMask{1},sliceAtlas{1},searchWin);
asd02 = asd(sliceMask{2},sliceAtlas{2},searchWin);
asd03 = asd(sliceMask{3},sliceAtlas{3},searchWin);

% Registered Mask and Atlas
dice11 = dice(logical(sliceTformedMask{1}),logical(sliceAtlas{1}));
dice12 = dice(logical(sliceTformedMask{2}),logical(sliceAtlas{2}));
dice13 = dice(logical(sliceTformedMask{3}),logical(sliceAtlas{3}));
hd11 = hd(sliceTformedMask{1},sliceAtlas{1},searchWin);
hd12 = hd(sliceTformedMask{2},sliceAtlas{2},searchWin);
hd13 = hd(sliceTformedMask{3},sliceAtlas{3},searchWin);
asd11 = asd(sliceTformedMask{1},sliceAtlas{1},searchWin);
asd12 = asd(sliceTformedMask{2},sliceAtlas{2},searchWin);
asd13 = asd(sliceTformedMask{3},sliceAtlas{3},searchWin);

% Dice
fprintf("Dice , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",dice01,dice02,dice03)
fprintf("Dice , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",dice11,dice12,dice13)

% Huasdorff
fprintf("\nHuasdorff , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",hd01,hd02,hd03)
fprintf("Huasdorff , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",hd11,hd12,hd13)

% Huasdorff
fprintf("\nASD , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",asd01,asd02,asd03)
fprintf("ASD , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",asd11,asd12,asd13)


% intersection of single Vertebras
for i = 1:4
    reg1 = pcReg{i};
    reg1 = reg1.Location;
    reg2 = pcReg{i+1};
    reg2 = reg2.Location;
    shp1 = alphaShape(reg1(:,1),reg1(:,2),reg1(:,3));
    shp2 = alphaShape(reg2(:,1),reg2(:,2),reg2(:,3));
    id1=inShape(shp2,reg1(:,1),reg1(:,2),reg1(:,3));
    id2=inShape(shp1,reg2(:,1),reg2(:,2),reg2(:,3));
    shp3=alphaShape([reg1(id1,1); reg2(id2,1)], [reg1(id1,2); reg2(id2,2)], ...
        [reg1(id1,3); reg2(id2,3)]);
    shpHolder{i} = shp3;
    V(i) = volume(shp3);
end
fprintf("\nVolumes:\n L1 & L2 = %.4f , L2 & L3 = %.4f , L3 & L4 = %.4f , L4 & L5 = %.4f\n",...
    V(1),V(2),V(3),V(4))

figure;
for i = 1:5
    reg = pcReg{i};
    reg = reg.Location;
    shp = alphaShape(reg(:,1),reg(:,2),reg(:,3));
    h = plot(shp);
    h.FaceAlpha = 0.05;
    hold on
end
for i = 1:4
    h = plot(shpHolder{i});
    h.FaceColor = 'r';
    hold on
end
title("Intersections After Single Vertebra CPD")

figure;
montage({sliceMask{1},sliceTformedMask{1},sliceAtlas{1},...
    sliceMask{2},sliceTformedMask{2},sliceAtlas{2},...
    sliceMask{3},sliceTformedMask{3},sliceAtlas{3}},'size',[3,3])
title("left: Mask - Middle: ICP Registered - Right: Atlas")



%% Single Vertebra Registeration - ICP then CPD
clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName = 'pat9'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


atlas = extractVertebra('pat0'); % can't touch this
atlasSave = atlas;
pat = extractVertebra(fileName);
patSave = pat;

% Extracting Outlayer
for i = 1:5
    [~,~,~,boundary] = outerLayer(atlas{i});
    atlas{i} = boundary;
    [~,~,~,boundary] = outerLayer(pat{i});
    pat{i} = boundary;
end

% ICP Registeration
for i = 1:5
    disp(i)
    [x,y,z] = ind2sub(size(atlas{i}),find(atlas{i}));
    pcTmpFixed = pointCloud([x,y,z]);
    pcTmpFixed.Normal = pcnormals(pcTmpFixed);
    [x,y,z] = ind2sub(size(pat{i}),find(pat{i}));
    pcTmpMoving = pointCloud([x,y,z]);
    pcTmpMoving.Normal = pcnormals(pcTmpMoving);
    [tformICP{i},pcRegICP{i}] = pcregistericp(pcTmpMoving,pcTmpFixed);
end

% CPD Registeration
for i = 1:5
    disp(i)
    [x,y,z] = ind2sub(size(atlas{i}),find(atlas{i}));
    pcTmpFixed = pointCloud([x,y,z]);
    pcTmpFixed = pcdownsample(pcTmpFixed,'nonuniformGridSample',12);
    pcTmpFixed.Count
    pcTmpMoving = pcdownsample(pcRegICP{i},'nonuniformGridSample',12);
    pcMovingBasePoints{i} = pcTmpMoving.Location;
    pcTmpMoving.Count
    [tform{i},pcReg{i}] = pcregistercpd(pcTmpMoving,pcTmpFixed,'verbose',true);
end

% ICP Interpolation for all points
pcRegSave = pcRegICP;
location = [];
for i = 1:5
    pcTmp = pcRegICP{i};
    location = [location ; pcTmp.Location];
end
pcMerged = pointCloud(location);

verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end
[xq,yq,zq] = ind2sub(size(mask),find(mask));
x = [];
y = [];
z = [];
for i = 1:5
    [xTmp,yTmp,zTmp] = ind2sub(size(pat{i}),find(pat{i}));
    x = [x ; xTmp];
    y = [y ; yTmp];
    z = [z ; zTmp];
end

regPoints = pcMerged.Location;
[x_tformed,y_tformed,z_tformed] = cpdInterpolation(x,y,z,regPoints,xq,yq,zq,'natural');

labelMatICP = ReconstructInterpMat(x_tformed,y_tformed,z_tformed,xq,yq,zq,mask);

% Removing NANs from ICP Interpolation
nanInds = find(isnan(x_tformed));
firstNanInds = nanInds;
x_tformed(nanInds) = [];
y_tformed(nanInds) = [];
z_tformed(nanInds) = [];

% CPD Interpolation
location = [];
for i = 1:5
    pcTmp = pcReg{i};
    location = [location ; pcTmp.Location];
end
pcMerged = pointCloud(location);
regPoints = pcMerged.Location;

xq = x_tformed;
yq = y_tformed;
zq = z_tformed;

location = [];
for i = 1:5
    pcTmp = pcMovingBasePoints{i};
    location = [location ; pcTmp];
end
x = location(:,1);
y = location(:,2);
z = location(:,3);

[x_tformed,y_tformed,z_tformed] = cpdInterpolation(x,y,z,regPoints,xq,yq,zq,'natural');

labelMat = ReconstructInterpMat(x_tformed,y_tformed,z_tformed,...
    round(xq),round(yq),round(zq),labelMatICP);


% Visulization
colorMat = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0; 0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880 ; 0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840];
% Atlas
figure;
subplot(1,3,2)
for i = 1:5
    [x,y,z] = ind2sub(size(atlasSave{i}),find(atlasSave{i}));
    shp = alphaShape(x,y,z);
    h = plot(shp);
    hold on
    h.FaceColor = colorMat(atlasSave{i}(x(1),y(1),z(1)),:);
end 
title("Atlas")

% Transformed
subplot(1,3,3)
for i = 1:5
    tmp = zeros(size(labelMat));
    tmp(find(labelMat == i+4)) = labelMat(find(labelMat == i+4));
    [x,y,z] = ind2sub(size(tmp),find(tmp));
    shp = alphaShape(x,y,z);
    h = plot(shp);
    h.FaceColor = colorMat(tmp(x(1),y(1),z(1)),:);
    hold on
end
title("Transformed Single Vertebras - CPD & ICP")

% Original
subplot(1,3,1)
for i = 1:5
    [x,y,z] = ind2sub(size(patSave{i}),find(patSave{i}));
    shp = alphaShape(x,y,z);
    h = plot(shp);
    h.FaceColor = colorMat(patSave{i}(x(1),y(1),z(1)),:);
    hold on
end 
title("Original Mask")



%% Registeration Quality Scores - Single Vertebra ICP The CPD
clc
tformPoints = [x_tformed,y_tformed,z_tformed];
nanInds = find(isnan(tformPoints(:,1)));

verteb = extractVertebra(fileName);
mask0 = 0;
for i = 1:5
    mask0 = mask0 + verteb{i};
end
[x0,y0,z0] = ind2sub(size(mask0),find(mask0));
x0(firstNanInds) = [];
y0(firstNanInds) = [];
z0(firstNanInds) = [];

% Jacobian
x_tformedJ = x_tformed;
y_tformedJ = y_tformed;
z_tformedJ = z_tformed;
nanIndsJ = find(isnan(x_tformed));
xq = x0;
yq = y0;
zq = z0;
xq(nanIndsJ) = [];
yq(nanIndsJ) = [];
zq(nanIndsJ) = [];
x_tformedJ(nanIndsJ) = [];
y_tformedJ(nanIndsJ) = [];
z_tformedJ(nanIndsJ) = [];
dx = x_tformedJ - xq;
dy = y_tformedJ - yq;
dz = z_tformedJ - zq;
sx = zeros(size(mask));
sy = zeros(size(mask));
sz = zeros(size(mask));
for i = 1:length(dx)
    sx(xq(i),yq(i),zq(i)) = dx(i);
    sy(xq(i),yq(i),zq(i)) = dy(i);
    sz(xq(i),yq(i),zq(i)) = dz(i);
end
jDetDisp(sx,sy,sz);


% Original Segmenation Mask
verteb = extractVertebra(fileName);
mask = 0;
for i = 1:5
    mask = mask + verteb{i};
end

% Atlas
verteb = extractVertebra('pat0');
atlas = 0;
for i = 1:5
    atlas = atlas + verteb{i};
end

% Transformed Segmentation Mask
tformedMask = labelMat;

% Dice Score
[~,sliceAtlas,indsAtlas] = extractSlice(atlas,atlas);
[~,sliceMask,indsMask] = extractSlice_ver2(mask,mask,indsAtlas);
[~,sliceTformedMask,indsTform] = extractSlice_ver2(tformedMask,tformedMask,indsAtlas);

searchWin = 43;
% Original Mask and Atlas
dice01 = dice(logical(sliceMask{1}),logical(sliceAtlas{1}));
dice02 = dice(logical(sliceMask{2}),logical(sliceAtlas{2}));
dice03 = dice(logical(sliceMask{3}),logical(sliceAtlas{3}));
hd01 = hd(sliceMask{1},sliceAtlas{1},searchWin);
hd02 = hd(sliceMask{2},sliceAtlas{2},searchWin);
hd03 = hd(sliceMask{3},sliceAtlas{3},searchWin);
asd01 = asd(sliceMask{1},sliceAtlas{1},searchWin);
asd02 = asd(sliceMask{2},sliceAtlas{2},searchWin);
asd03 = asd(sliceMask{3},sliceAtlas{3},searchWin);

% Registered Mask and Atlas
dice11 = dice(logical(sliceTformedMask{1}),logical(sliceAtlas{1}));
dice12 = dice(logical(sliceTformedMask{2}),logical(sliceAtlas{2}));
dice13 = dice(logical(sliceTformedMask{3}),logical(sliceAtlas{3}));
hd11 = hd(sliceTformedMask{1},sliceAtlas{1},searchWin);
hd12 = hd(sliceTformedMask{2},sliceAtlas{2},searchWin);
hd13 = hd(sliceTformedMask{3},sliceAtlas{3},searchWin);
asd11 = asd(sliceTformedMask{1},sliceAtlas{1},searchWin);
asd12 = asd(sliceTformedMask{2},sliceAtlas{2},searchWin);
asd13 = asd(sliceTformedMask{3},sliceAtlas{3},searchWin);

% Dice
fprintf("Dice , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",dice01,dice02,dice03)
fprintf("Dice , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",dice11,dice12,dice13)

% Huasdorff
fprintf("\nHuasdorff , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",hd01,hd02,hd03)
fprintf("Huasdorff , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",hd11,hd12,hd13)

% Huasdorff
fprintf("\nASD , Original Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",asd01,asd02,asd03)
fprintf("ASD , Registered Mask & Atlas:\n Lateral: %.4f Axial: %.4f Frontal: %.4f \n",asd11,asd12,asd13)


% intersection of single Vertebras
for i = 1:4
    reg1 = pcReg{i};
    reg1 = reg1.Location;
    reg2 = pcReg{i+1};
    reg2 = reg2.Location;
    shp1 = alphaShape(reg1(:,1),reg1(:,2),reg1(:,3));
    shp2 = alphaShape(reg2(:,1),reg2(:,2),reg2(:,3));
    id1=inShape(shp2,reg1(:,1),reg1(:,2),reg1(:,3));
    id2=inShape(shp1,reg2(:,1),reg2(:,2),reg2(:,3));
    shp3=alphaShape([reg1(id1,1); reg2(id2,1)], [reg1(id1,2); reg2(id2,2)], ...
        [reg1(id1,3); reg2(id2,3)]);
    shpHolder{i} = shp3;
    V(i) = volume(shp3);
end
fprintf("\nVolumes:\n L1 & L2 = %.4f , L2 & L3 = %.4f , L3 & L4 = %.4f , L4 & L5 = %.4f\n",...
    V(1),V(2),V(3),V(4))

figure;
for i = 1:5
    reg = pcReg{i};
    reg = reg.Location;
    shp = alphaShape(reg(:,1),reg(:,2),reg(:,3));
    h = plot(shp);
    h.FaceAlpha = 0.05;
    hold on
end
for i = 1:4
    h = plot(shpHolder{i});
    h.FaceColor = 'r';
    hold on
end
title("Intersections After Single Vertebra CPD")

figure;
montage({sliceMask{1},sliceTformedMask{1},sliceAtlas{1},...
    sliceMask{2},sliceTformedMask{2},sliceAtlas{2},...
    sliceMask{3},sliceTformedMask{3},sliceAtlas{3}},'size',[3,3])
title("left: Mask - Middle: ICP Registered - Right: Atlas")




