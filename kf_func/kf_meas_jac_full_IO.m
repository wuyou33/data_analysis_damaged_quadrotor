function Hx = kf_meas_jac_full_IO(x,u,dt)


% Hx = [diag(ones(6,1)) zeros(6,9)];

p = u(1) - x(7);
q = u(2) - x(8);
r = u(3) - x(9);

A = [0 -r q; r 0 -p; -q p 0];
Hx_att = [diag(ones(3,1)) zeros(3,12)];
Hx_v   = [zeros(3,3) diag(ones(3,1)) zeros(3,6) A];

Hx = [Hx_att;Hx_v];
return