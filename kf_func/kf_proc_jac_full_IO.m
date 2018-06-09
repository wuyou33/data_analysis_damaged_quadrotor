function Fx = kf_proc_jac_full_IO(x,u,dt)

p = u(1) - x(7);
q = u(2) - x(8);
r = u(3) - x(9);
ax = u(4) - x(10);
ay = u(5) - x(11);
az = u(6) - x(12);
phi = x(1);
theta = x(2);
%psi = x(3);
vx = x(4);
vy = x(5);
vz = x(6);
p_dot = u(7);
q_dot = u(8);
r_dot = u(9);
%r1 r2 r3 indicate displacement from cg to IMU

g = 9.81;

Fx = zeros(15);

%x_k(1) = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
Fx(1,1) = q*cos(phi)*tan(theta) - r*sin(phi)*tan(theta);
Fx(1,2) = q*sin(phi)*(sec(theta))^2 + r*cos(phi)*(sec(theta))^2;
Fx(1,7) = -1;
% added later:
Fx(1,8) = -(sin(phi)*tan(theta));
Fx(1,9) = -(cos(phi)*tan(theta));

% Fx(2,1) = -r*cos(phi);
% Fx(2,2) = -q*sin(theta);
% Fx(2,8) = -cos(theta);
% Fx(2,9) = sin(phi);
% CORRECTION
%x_k(2) = q*cos(phi) - r*sin(phi); % correct equation!
Fx(2,1) = -q*sin(phi) - r*cos(phi);
Fx(2,8) = -cos(phi);
Fx(2,9) = sin(phi);

%x_k(3) = q*sin(phi)*sec(theta) + r*cos(phi)*sec(theta);
Fx(3,1) = q*cos(phi)*sec(theta) - r*sin(phi)*sec(theta);
Fx(3,2) = q*sin(phi)*sec(theta)*tan(theta) + r*cos(phi)*sec(theta)*tan(theta);
Fx(3,8) = -sin(phi)*sec(theta);
Fx(3,9) = -cos(phi)*sec(theta);

%x_k(4) = ax - q*vz + r*vy - g*sin(theta);
Fx(4,2) = -g*cos(theta);
Fx(4,5) = r;
Fx(4,6) = -q;
Fx(4,8) = vz;
Fx(4,9) = -vy;
Fx(4,10) = -1;
Fx(4,13) = q^2+r^2;
Fx(4,14) = -p*q+r_dot;
Fx(4,15) = -p*r-q_dot;

%x_k(5) = ay - r*vx + p*vz + g*sin(phi)*cos(theta);
Fx(5,1) = g*cos(phi)*cos(theta);
Fx(5,2) = -g*sin(phi)*sin(theta);
Fx(5,4) = -r;
Fx(5,6) = p; 
Fx(5,7) = -vz;
Fx(5,9) = vx;
Fx(5,11) = -1;
Fx(5,13) = -p*q-r_dot;
Fx(5,14) = p^2+r^2;
Fx(5,15) = -q*r+p_dot;

%x_k(6) = az - p*vy + q*vx + g*cos(theta)*cos(phi);
Fx(6,1) = -g*cos(theta)*sin(phi);
Fx(6,2) = -g*sin(theta)*cos(phi);
Fx(6,4) = q;
Fx(6,5) = -p;
Fx(6,7) = vy;
Fx(6,8) = -vx;
Fx(6,12) = -1;
Fx(6,13) = -p*r+q_dot;
Fx(6,14) = -q*r-p_dot;
Fx(6,15) = p^2+q^2;

% x_k(1) = ax - q*x(3) + r*x(2) - g*sin(theta);
% x_k(2) = ay - r*x(1) + p*x(3) + g*sin(phi)*cos(theta);
% x_k(3) = az - p*x(2) + q*x(1) + g*cos(theta)*cos(phi);


return
