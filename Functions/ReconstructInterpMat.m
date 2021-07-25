function labelMat = ReconstructInterpMat(x_tformed,y_tformed,z_tformed,xq,yq,zq,vertebMask)

% getting rid of NANs
nanInds = find(isnan(x_tformed));
xTmp = x_tformed;
xTmp(nanInds) = [];
yTmp = y_tformed;
yTmp(nanInds) = [];
zTmp = z_tformed;
zTmp(nanInds) = [];

xqTmp = xq;
xqTmp(nanInds) = [];
yqTmp = yq;
yqTmp(nanInds) = [];
zqTmp = zq;
zqTmp(nanInds) = [];

% Reconstruct Label Matrix
sz = size(vertebMask);
labelMat = zeros(sz);
for i = 1:length(xTmp)
    ind1 = round(xTmp(i));
    ind2 = round(yTmp(i));
    ind3 = round(zTmp(i));
    labelMat(ind1,ind2,ind3) = vertebMask(xqTmp(i),yqTmp(i),zqTmp(i));
end
