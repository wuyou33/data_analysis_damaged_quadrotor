% Rotation around x-axes

function S=Rot_x(fix);

S=[1 0        0
   0 cos(fix) -sin(fix)
   0 sin(fix) cos(fix)];