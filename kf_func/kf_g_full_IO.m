function G = kf_g_full_IO(x)
phi = x(1);
theta = x(2);
%psi = x(3);
vx = x(4);
vy = x(5);
vz = x(6);

g = 9.81;

G = zeros(15,6);

G(1,1) = -1;
G(1,2) = -sin(phi)*tan(theta);
G(1,3) = -cos(phi)*tan(theta);

G(2,2) = -cos(phi);
G(2,3) = sin(phi);

G(3,2) = -sin(phi)*sec(theta);
G(3,3) = -cos(phi)*sec(theta);

G(4,2) = vz; 
G(4,3) = -vy;
G(4,4) = -1; 

G(5,1) = -vz; 
G(5,3) = vx; 
G(5,5) = -1;

G(6,1) = vy;
G(6,2) = -vx;
G(6,6) = -1;

 return