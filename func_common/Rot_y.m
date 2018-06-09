% Rotation around y-axes

function S=Rot_y(fiy);

S=[cos(fiy)  0 sin(fiy)
   0         1 0       
   -sin(fiy) 0 cos(fiy)];