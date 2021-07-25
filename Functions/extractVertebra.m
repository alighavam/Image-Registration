function verteb = extractVertebra(fileName)

[~,mask,sz] = loadData(fileName);
mask = flip(mask,3);
regLabels = [20,21,22,23,24]-15;
for i = 1:length(regLabels)
    label = regLabels(i);
    tmp = zeros(sz);
    tmp(find(mask == label)) = mask(find(mask == label));
    verteb{i} = tmp;
end