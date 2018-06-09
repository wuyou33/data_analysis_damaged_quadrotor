% Rotation around z-axes

function S=Rot_z(fiz);

S=[cos(fiz) -sin(fiz) 0
   sin(fiz) cos(fiz)  0
   0        0         1];