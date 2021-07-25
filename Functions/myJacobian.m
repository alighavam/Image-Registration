function det_J = myJacobian(sx,sy,sz)

[gx_y,gx_x,gx_z] = gradient(sx);
[gy_y,gy_x,gy_z] = gradient(sy);
[gz_y,gz_x,gz_z] = gradient(sz);

gx_x = gx_x + 1;
gy_y = gy_y + 1;
gz_z = gz_z + 1;

det_J = gx_x.*gy_y.*gz_z + ...
        gy_x.*gz_y.*gx_z + ...
        gz_x.*gx_y.*gy_z - ...
        gz_x.*gy_y.*gx_z - ...
        gy_x.*gx_y.*gz_z - ...
        gx_x.*gz_y.*gy_z;

% Source
% https://www.programmersought.com/article/89892588427/