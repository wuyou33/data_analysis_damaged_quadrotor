psi = 0:1:720;

w1obs = interp1(psi_mod2(du)-52.7,OB_a.w1obs(du),psi,'near');
w2obs = interp1(psi_mod2(du)+52.7,OB_a.w2obs(du),psi,'near');
w3obs = interp1(psi_mod2(du)+127.3,OB_a.w3obs(du),psi,'near');
w4obs = interp1(psi_mod2(du)-127.3,OB_a.w4obs(du),psi,'near');

w1ref = interp1(psi_mod2(du)-52.7,OB_a.w1ref(du),psi,'near');
w2ref = interp1(psi_mod2(du)+52.7,OB_a.w2ref(du),psi,'near');
w3ref = interp1(psi_mod2(du)+127.3,OB_a.w3ref(du),psi,'near');
w4ref = interp1(psi_mod2(du)-127.3,OB_a.w4ref(du),psi,'near');

r  = interp1(psi_mod2(du),OB_a.R(du),psi,'near');

x_grid = cos((psi)/57.3);
y_grid = sin((psi)/57.3);


% figure
% plot3(y_grid,x_grid,w1obs,'.'); hold on;
% plot3(y_grid,x_grid,w2obs,'.');
% plot3(y_grid,x_grid,w4obs,'.');
% xlabel('y'); ylabel('x'); grid on;

figure
plot3(y_grid,x_grid,w1ref.^2+w2ref.^2+w4ref.^2,'b.'); hold on;
plot3(y_grid,x_grid,w1obs.^2+w2obs.^2+w4obs.^2,'r.'); hold on;
xlabel('y'); ylabel('x'); grid on;title('sum of \omega'); legend('ref','obs')

% figure
% plot3(y_grid,x_grid)