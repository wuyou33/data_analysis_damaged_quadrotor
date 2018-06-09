function x_k = kf_proc_full_IO(x,u,dt)

%r_u1 = x(13);
%r_u2 = x(14);
%r_u3 = x(15);

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

g = 9.81;

x_k = zeros(1,12);

x_k(1) = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
x_k(2) = q*cos(phi) - r*sin(phi);
x_k(3) = q*sin(phi)*sec(theta) + r*cos(phi)*sec(theta);
x_k(4) = ax - q*vz + r*vy - g*sin(theta);
x_k(5) = ay - r*vx + p*vz + g*sin(phi)*cos(theta);
x_k(6) = az - p*vy + q*vx + g*cos(theta)*cos(phi);
x_k(7) = 0;
x_k(8) = 0;
x_k(9) = 0;
x_k(10) = 0;
x_k(11) = 0;
x_k(12) = 0;

x_k(4) = x_k(4) + (q^2+r^2)*x(13) - p*q*x(14) - p*r*x(15) + r_dot*x(14) - q_dot*x(15);
x_k(5) = x_k(5) - p*q*x(13) + (p^2+r^2)*x(14) - q*r*x(15) - r_dot*x(13) + p_dot*x(15);
x_k(6) = x_k(6) - p*r*x(13) - q*r*x(14) + (p^2+q^2)*x(15) + q_dot*x(13) - p_dot*x(14);

x_k(13) = 0;
x_k(14) = 0;
x_k(15) = 0;
return
