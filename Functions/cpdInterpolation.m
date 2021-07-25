function [x_tformed,y_tformed,z_tformed] = cpdInterpolation(x,y,z,regPoints,xq,yq,zq,method)

% interpolation for X dimension --> x_tformed = f1(x,y,z), we will find f1
tic
x_tformed = griddata(x,y,z,regPoints(:,1),xq,yq,zq,method);
toc

% interpolation for Y dimension --> y_tformed = f2(x,y,z), we will find f2
tic
y_tformed = griddata(x,y,z,regPoints(:,2),xq,yq,zq,method);
toc

% interpolation for Z dimension --> z_tformed = f3(x,y,z), we will find f3
tic
z_tformed = griddata(x,y,z,regPoints(:,3),xq,yq,zq,method);
toc