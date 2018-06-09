function y_k = kf_meas_full_IO(x,u,dt)

% y_k = (x(1:6))';
%rw indicates diplacement from Optitrack origin to IMU. User defined.
rw = [0.003 -0.005 0]; 

p = u(1) - x(7);
q = u(2) - x(8);
r = u(3) - x(9);
vx = x(4);
vy = x(5);
vz = x(6);
r1 = x(13)-rw(1);
r2 = x(14)-rw(2);
r3 = x(15)-rw(3);

y_k(1) = x(1);
y_k(2) = x(2);
y_k(3) = x(3);
y_k(4) = vx - r*r2 + q*r3;
y_k(5) = vy + r*r1 - p*r3;
y_k(6) = vz - q*r1 + p*r2;

return
