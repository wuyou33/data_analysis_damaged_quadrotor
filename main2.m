
clear all

for  index = 47; %b47 r54 81

addpath(genpath('E:\Data\_code_import_files'));
[OB_a,OT_a,WIND,PARA,take,DU] = import_data(index,1, 1, 1,1,1,'Q',0);

%% plot raw data
% figure
% subplot(2,1,1)
% plot(OB_a.TIME,OB_a.h1); hold on; grid on; title(['h1,h2',' ',take.name]);
% % xlim([22,24]);
% subplot(2,1,2)
% plot(OB_a.TIME,OB_a.h2); grid on; 
% 
% w1ref_filter = butterworth(OB_a.w1ref,4,50/256);
% w2ref_filter = butterworth(OB_a.w2ref,4,50/256);
% w3ref_filter = butterworth(OB_a.w3ref,4,50/256);
% w4ref_filter = butterworth(OB_a.w4ref,4,50/256);
% figure; set(gcf,'position',[0,0,400,800]);
% subplot(4,1,1)
% plot(OB_a.TIME,w1ref_filter); hold on; grid on; title(['rotor speed',' ',take.name]);
% plot(OB_a.TIME,OB_a.w1obs);
% % xlim([22,24]);
% subplot(4,1,2)
% plot(OB_a.TIME,w2ref_filter); hold on; grid on; title(['rotor speed',' ',take.name]);
% plot(OB_a.TIME,OB_a.w2obs);
% % xlim([22,24]);
% subplot(4,1,3)
% plot(OB_a.TIME,w3ref_filter); hold on; grid on; title(['rotor speed',' ',take.name]);
% plot(OB_a.TIME,OB_a.w3obs);
% % xlim([22,24]);
% subplot(4,1,4)
% plot(OB_a.TIME,w4ref_filter); hold on; grid on; title(['rotor speed',' ',take.name]);
% plot(OB_a.TIME,OB_a.w4obs);
% % xlim([22,24]);
% 
% figure
% subplot(3,1,1)
% plot(OB_a.TIME,OB_a.p_des); hold on; grid on; title(['pq_{des} and pq',' ',take.name]);
% plot(OB_a.TIME,OB_a.p);
% % xlim([22,24]);
% % plot(OT_a.TIME,OT_a.P);
% subplot(3,1,2)
% plot(OB_a.TIME,OB_a.q_des); hold on; grid on;
% plot(OB_a.TIME,OB_a.q);
% % xlim([22,24]);
% % plot(OT_a.TIME,OT_a.Q);
% subplot(3,1,3)
% plot(OB_a.TIME,OB_a.r,'r'); hold on; grid on;
%% preprocessing data

DU = [1, length(OT_a.TIME)];
DU = [round(length(OT_a.TIME)/20), round(length(OT_a.TIME)*9.5/10)];
DU = [10920 139600]; %RB n = 0;
% DU = [1 137900]; %LB n = 0
% DU = [24950 32450]; %#48 v = 1
% DU = [47680 51400]; %#48 v = 2
% DU = [62510 67440]; %#48 v = 3
% DU = [79030 82970]; %#48 v = 4
% DU = [91720 95970]; %#48 v = 5
% DU = [105700 109900]; %#48 v = 6
% DU = [120300 123900]; %#48 v = 7
% DU = [132600 139300]; %#48 v = 8

% DU = [9729 13260]; %#47 v = 1
% DU = [28250 32120]; %#47 v = 2
% DU = [43380 46690]; %#47 v = 3
DU = [57080 60620]; %#47 v = 4
% DU = [73450 76260]; %#47 v = 5
% DU = [83860 87020]; %#47 v = 6
% DU = [93300 96690]; %#47 v = 7
% DU = [105400 106000]; %#47 v = 8

% DU = [7240 12000]; %#50 v = 1
% DU = [26970 31800]; %#50 v = 2
% DU = [41920 47840]; %#50 v = 3
% DU = [60170 64050]; %#50 v = 4
% DU = [80970 76550]; %#50 v = 5
% DU = [90370 94510]; %#50 v = 6
% DU = [106400 110900]; %#50 v = 7
% DU = [121600 124600]; %#50 v = 8
% DU = [130500 131500]; %#50 v = 9

fc = 3.0;
du = DU(1):DU(2);
% du = 76790:82910;
time = (0:length(du)-1)/512;

% calculate acceleartion from external sensor (Optitrack)
VX_filt = butterworth(OT_a.VX,4,fc/256);
VY_filt = butterworth(OT_a.VY,4,fc/256);
VZ_filt = butterworth(OT_a.VZ,4,fc/256);
VX_filt = butterworth(VX_filt,4,fc/256);
VY_filt = butterworth(VY_filt,4,fc/256);
VZ_filt = butterworth(VZ_filt,4,fc/256);




wa= (OB_a.w1obs+OB_a.w2obs+OB_a.w3obs+OB_a.w4obs)/4;
VY_air_filt = butterworth(OT_a.VY_air,4,fc/256);
Vz_filt = butterworth(OT_a.VZ,4,4/256);
azi_filt = derivative(Vz_filt,OB_a.TIME);
psi_mod = mod(OT_a.PSI+90,360);
psi_mod2 = mod(OT_a.PSI,720);

R = 0.075;
w_bar = sqrt((OB_a.w1obs.^2+OB_a.w2obs.^2+OB_a.w3obs.^2)/3)*2*pi/60;
Area = 3*pi*R^2;
rho = 1.225;

V_filt = sqrt(VX_filt.^2+VY_air_filt.^2+VZ_filt.^2);

% du1 = find(VY_air_filt<=-3.90);
% du2 = find(VY_air_filt>=-4.10);
% du = intersect(du1,du2);
%% analyze states in different heading angle

r = compress_matrix(OB_a.R(du),50);
w1 = compress_matrix(OB_a.w1obs(du),50)*2*pi/60;
w2 = compress_matrix(OB_a.w2obs(du),50)*2*pi/60;
w3 = compress_matrix(OB_a.w3obs(du),50)*2*pi/60;
psi_mod_comp = compress_matrix(psi_mod(du),50);


% figure(4)
% plot(psi_mod(du),r,'.'); title(take.name); hold on;
% plot(psi_mod_comp,r,'.'); title(take.name); hold on;

figure(11)
subplot(2,3,5)
plot(psi_mod_comp,w1,'+'); hold on;
plot(psi_mod_comp,w2,'o');
plot(psi_mod_comp,w3,'x');
subplot(1,3,3)
plot(psi_mod_comp,r,'+'); hold on;
% continue;
%% Stability frame calculation. Force on the stability frame.


dVX = derivative(VX_filt,OT_a.TIME);
dVY = derivative(VY_filt,OT_a.TIME);
dVZ = derivative(VZ_filt,OT_a.TIME);

% project n on Inertia frame; project velocity and dV on Body frame. 
na_I = zeros(length(OT_a.TIME),2); 

n = [0.0 0.0 -1.0]';
% n   = [0.155 -0.005 -1.0]'; %primary axis

n_I = zeros(length(OT_a.TIME),3);
Ab  = zeros(length(OT_a.TIME),3); 
Vb  = zeros(length(OT_a.TIME),3); 
for i = 1:length(OT_a.TIME)
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
    
    R_IB = R_BI';
    
    Ab(i,:) = [dVX(i) dVY(i) dVZ(i)]*R_BI';
    Vb(i,:) = [VX_filt(i) VY_filt(i) VZ_filt(i)]*R_BI';
    n_I(i,:) =  n'*R_IB';
end

% calculate Stability frame. Project Velocity and dV on Stability frame.
na_I(:,1) = butterworth(n_I(:,1),4,fc/512);
na_I(:,2) = butterworth(n_I(:,2),4,fc/512);
% na_I(:,1) = n_I(:,1);
% na_I(:,2) = n_I(:,2);

phi_s = asin(na_I(:,2));
theta_s = -asin(na_I(:,1)./cos(phi_s));


As = zeros(length(OT_a.TIME),3);
Vs = zeros(length(OT_a.TIME),3);
for i = 1:length(OT_a.TIME)
    theta = theta_s(i); phi = phi_s(i);
    R_SI = [cos(theta) 0 -sin(theta);
            sin(phi)*sin(theta) cos(phi) sin(phi)*cos(theta);
            cos(phi)*sin(theta) -sin(phi) cos(phi)*cos(theta)];    
        
    As(i,:) = [dVX(i) dVY(i) dVZ(i)-9.81]*R_SI';
    Vs(i,:) = [VX_filt(i) VY_air_filt(i) VZ_filt(i)] * R_SI';    
end


Vs_h = sqrt(Vs(:,1).^2+Vs(:,2).^2);
As_h = sqrt(As(:,1).^2+As(:,2).^2);


beta_v = asin(Vs(:,1)./Vs_h);
beta_a = asin(As(:,1)./As_h);

psi_mod = mod(OT_a.PSI,360)-180;

% plot figures

% figure
% plot(OT_a.TIME,phi_s);
% hold on
% plot(OT_a.TIME,theta_s);
% xlabel('TIME'); legend('phi_s','theta_s');
% 
% figure
% subplot(3,1,1);
% plot(Vs(du,1),As(du,1));  ylabel('As_x');title(take.name);
% grid on;
% subplot(3,1,2);
% plot(Vs(du,2),As(du,2)/0.81);grid on; ylabel('As_y');
% subplot(3,1,3);
% plot(Vs(du,3),As(du,3));grid on;xlabel('Vy_{air}'); ylabel('As_z');
% 

figure
subplot(3,1,1);
plot(sqrt(OT_a.VX_air(du).^2+VY_air_filt(du).^2+OT_a.VZ_air(du).^2),As(du,1));  ylabel('As_x');title(take.name);
grid on;
subplot(3,1,2);
plot(VY_air_filt(du),As(du,2));grid on; ylabel('As_y');
subplot(3,1,3);
plot(VY_air_filt(du),As(du,3));grid on;xlabel('Vy_{air}'); ylabel('As_z');

% figure
% subplot(3,1,2);
% plot(-VY_air_filt(du),-As(du,1));  ylabel('As_y');
% grid on;
% subplot(3,1,1);
% plot(-VY_air_filt(du),As(du,2));grid on; ylabel('As_x');title(take.name);
% subplot(3,1,3);
% plot(-VY_air_filt(du),As(du,3));grid on;xlabel('Vy_{air}'); ylabel('As_z');

% figure
% plot(VY_air_filt(du),OB_a.R(du)); ylabel('R'); xlabel('Vy_{air}');title(take.name);
% % 
% 
% figure
% plot(psi_mod(du),OB_a.R(du));
% 
% figure
% plot(Vs(du,2),As(du,2)*PARA.mass/1000); grid on; xlabel('Vs_y'); ylabel('As_y');


% figure
% plot(Vs(du,2),As(du,3)*PARA.mass/1000); grid on; xlabel('Vs_y'); ylabel('Fz');

% figure
% plot(Vs(du,2),-As(du,3)*PARA.mass/1000./w_bar(du).^2/Area/R^2/rho); grid on; xlabel('Vs_y'); ylabel('Ct')

%% Aerodynamc moment, body axis.
Ip = PARA.Ip;
Iv = PARA.Iv;

k0 = 1.9035e-6;
t0 = 1.9202951e-8;
dr = -1.918988e-3;

         
w1 = OB_a.w1obs(du)*2*pi/60;
w2 = OB_a.w2obs(du)*2*pi/60;
w3 = OB_a.w3obs(du)*2*pi/60;
w4 = 0*OB_a.w4ref(du)*2*pi/60;

[b,a] = butter(2,16/256);
p =  filtfilt(b,a,OB_a.p(du));
q =  filtfilt(b,a,OB_a.q(du));
r =  filtfilt(b,a,OB_a.r(du));
p_dot = [0; diff(p)]*512;
q_dot = [0; diff(q)]*512;
r_dot = [0; diff(r)]*512;
w1_dot = [0; diff(w1)]*512;
w2_dot = [0; diff(w2)]*512;
w3_dot = [0; diff(w3)]*512;

Mx = p_dot *Iv(1,1) - (Iv(2,2)-Iv(3,3)).*q .*r  ...
    + Ip(3,3).*q .*w1  - Ip(3,3).*q .*w2 ...
    + Ip(3,3).*q .*w3 - Ip(3,3).*q .*w4;

My = q_dot *Iv(2,2) - (Iv(3,3)-Iv(1,1)).*p .*r  ...
    - Ip(3,3).*p .*w1 + Ip(3,3).*p .*w2 ...
    - Ip(3,3).*p .*w3 + Ip(3,3).*p .*w4;

Mz = r_dot *Iv(3,3) - (Iv(1,1)-Iv(2,2)).*p .*q  ...
     - w1_dot *Ip(3,3) + w2_dot *Ip(3,3) - w3_dot *Ip(3,3);

up = w1.^2-w2.^2-w3.^2;
uq = w1.^2+w2.^2-w3.^2;
ur = w1.^2-w2.^2+w3.^2;

l = 0.0875;
b = 0.115;

dMx = Mx - k0*b*up;
dMy = My - k0*l*uq;
dMz = Mz - t0*ur - dr*r;

figure(11)
subplot(2,3,5)
plot(compress_matrix(psi_mod(du),50),compress_matrix(dMx,50),'+'); hold on;
plot(compress_matrix(psi_mod(du),50),compress_matrix(dMy,50),'o'); hold on;

% figure
% plot3(compress_matrix(OT_a.VY_air,20),compress_matrix(psi_mod,20),compress_matrix(dMx,20),'.')
end
return;


%% Plot

% figure
% subplot(2,1,1)
% plot(OB_a.TIME,ax_filt); hold on;
% plot(OB_a.TIME,ay_filt);
% subplot(2,1,2)
% plot(OB_a.TIME,Ax_g); hold on;
% plot(OB_a.TIME,Ay_g);
% 
% figure
% subplot(2,3,1);
% pwelch(OB_a.ax,[],[],[],512);
% subplot(2,3,2);
% pwelch(OB_a.ay,[],[],[],512);
% subplot(2,3,3);
% pwelch(OB_a.az,[],[],[],512);
% subplot(2,3,4);
% pwelch(OB_a.p,[],[],[],512);
% subplot(2,3,5);
% pwelch(OB_a.q,[],[],[],512);
% subplot(2,3,6);
% pwelch(OB_a.r,[],[],[],512);


% xlim([22,24]);
% fs = 512;
% figure
% subplot(2,1,1)
% periodogram(OB_a.ax,[],[],fs);
% subplot(2,1,2)
% periodogram(OB_a.ay,[],[],fs);
% 
% return;

figure
subplot(3,1,1)
plot(OB_a.TIME,OB_a.phi*57.3); hold on; grid on;title(['attitude',' ',take.name]);
% plot(OT_a.TIME-OT_a.TIME(1),OT_a.PHI);
plot(OB_a.TIME,OB_a.phi_ot*57.3);
% xlim([22,24]);
subplot(3,1,2)
plot(OB_a.TIME,OB_a.theta*57.3); hold on; grid on
% plot(OT_a.TIME-OT_a.TIME(1),OT_a.THETA);
plot(OB_a.TIME,OB_a.theta_ot*57.3);
% xlim([22,24]);
subplot(3,1,3)
plot(OB_a.TIME,OB_a.psi*57.3); hold on; grid on
% plot(OT_a.TIME-OT_a.TIME(1),mod(OT_a.PSI,360)-180);
plot(OB_a.TIME,OB_a.psi_ot*57.3);
% xlim([22,24]);


% figure
% subplot(3,1,1)
% plot(OB_a.TIME,OB_a.ax);hold on; grid on;title(['az^b',' ',take.name]);
% plot(OB_a.TIME,ax_filt);
% subplot(3,1,2)
% plot(OB_a.TIME,OB_a.ay);hold on; grid on;
% plot(OB_a.TIME,ay_filt);
% subplot(3,1,3)
% plot(OB_a.TIME,OB_a.az);hold on; grid on;
% plot(OB_a.TIME,az_filt);
% 

figure
% subplot(3,1,1)
% plot(OB_a.TIME,OB_a.p);hold on; grid on;title('pqr');
% subplot(3,1,2)
% plot(OB_a.TIME,OB_a.q);hold on; grid on;
% subplot(3,1,3)
% plot(OB_a.TIME,OB_a.r);hold on; grid on;



% figure
% plot(OB_a.TIME,OB_a.w1obs_indi); hold on; grid on;title(['w',' ',take.name]);
% plot(OB_a.TIME,OB_a.w2obs_indi);
% plot(OB_a.TIME,OB_a.w3obs_indi);
% plot(OB_a.TIME,OB_a.w4obs_indi);

figure
subplot(3,1,1)
plot(OB_a.TIME,OB_a.w1ref); hold on; grid on; title(['rotor speed',' ',take.name]);
plot(OB_a.TIME,OB_a.w1obs);
% xlim([22,24]);
subplot(3,1,2)
plot(OB_a.TIME,OB_a.w2ref); hold on; grid on; title(['rotor speed',' ',take.name]);
plot(OB_a.TIME,OB_a.w2obs);
% xlim([22,24]);
subplot(3,1,3)
plot(OB_a.TIME,OB_a.w3ref); hold on; grid on; title(['rotor speed',' ',take.name]);
plot(OB_a.TIME,OB_a.w3obs);
% xlim([22,24]);


%%
% 
figure
plot(OB_a.TIME,OB_a.r); hold on; grid on; title('r comparison');
plot(OB_a.TIME,OB_a.r_ot);
% plot(OB_a.TIME,OB_a.psi_ot*10);
% plot(OT_a.TIME,OT_a.PSI);

% figure
% subplot(2,1,1)
% plot(OB_a.TIME,sqrt(OT_a.U.^2+OT_a.V.^2));
% ylabel('V_h (m/s)')
% subplot(2,1,2)
% plot(OB_a.TIME,sqrt(OT_a.U.^2+OT_a.V.^2+OT_a.W.^2));
% ylabel('V (m/s)')
% 
figure
subplot(3,1,1)
plot(OB_a.TIME,OB_a.acc_des_x); hold on;
plot(OB_a.TIME,OB_a.acc_des_x_filter); grid on;
subplot(3,1,2)
plot(OB_a.TIME,OB_a.acc_des_y); hold on;
plot(OB_a.TIME,OB_a.acc_des_y_filter); grid on;
subplot(3,1,3)
plot(OB_a.TIME,OB_a.acc_des_z); hold on;
plot(OB_a.TIME,OB_a.acc_des_z_filter); grid on;

% figure
% subplot(2,1,1)
% plot(OB_a.TIME, OB_a.p_des_dot); hold on; grid on;
% subplot(2,1,2)
% plot(OB_a.TIME, OB_a.q_des_dot); hold on; grid on;
% % 
% figure
% subplot(2,1,1)
% plot(OB_a.TIME, OB_a.p_des); hold on; grid on;
% plot(OB_a.TIME, OB_a.p_des_filter); 
% subplot(2,1,2)
% plot(OB_a.TIME, OB_a.q_des); hold on; grid on;
% plot(OB_a.TIME, OB_a.q_des_filter); 

figure
plot3(-OT_a.X,OT_a.Y,-OT_a.Z); grid on; xlabel('X'); ylabel('Y'); zlabel('H');

figure
plot(OB_a.TIME,OT_a.VX); grid on; hold on;
plot(OB_a.TIME,OT_a.VY); grid on; legend('Vx','Vy')
% return;

%% Relationship between r and n tilt
h1 = OB_a.h1; h2 = OB_a.h2;
hh = sqrt(abs(h1.^2+h2.^2));
ksi = asin(hh);
%% Analysis Vs vs Ct
W1 = OB_a.w1obs; W2 = OB_a.w2obs; W3 = OB_a.w3obs; W4 = 0*OB_a.w4obs;
W_bar = sqrt((W1.^2+W2.^2+W3.^2+W4.^2)/4);
w1 = W1./W_bar; w2 = W2./W_bar; w3 = W3./W_bar; w4 = W4./W_bar;
R = 0.06;
Area = pi.*R.^2;
rho = 1.225;
m = 0.5;

T = az_filt*m;
Ct = T(du)./W_bar(du).^2/Area/R^2/rho;

W_bar_filt = butterworth(W_bar,4,fc/256);
T_filt = As(:,3)*m;
Ct_filt = T(du)./W_bar_filt(du).^2/Area/R^2/rho;

%% solve rw by differentiating Vo = Vc + w x ru
VX_filt = butterworth(OT_a.VX,4,10/512);
VY_filt = butterworth(OT_a.VY,4,10/512);
VZ_filt = butterworth(OT_a.VZ,4,10/512);
dU = derivative(VX_filt,OT_a.TIME);
dV = derivative(VY_filt,OT_a.TIME);
dW = derivative(VZ_filt,OT_a.TIME);
dP = derivative(OB_a.P,OT_a.TIME);
dQ = derivative(OB_a.Q,OT_a.TIME);
dR = derivative(OB_a.R,OT_a.TIME);

rw = [0.003 -0.005 0];
% rw = zeros(3,length(OT_a.TIME));
da = zeros(3,length(OT_a.TIME));
datest = zeros(3,length(OT_a.TIME)); 
dVo_B = zeros(3,length(OT_a.TIME)); 
g_B = zeros(3,length(OT_a.TIME)); 
for i = 1:length(OT_a.TIME)
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
    dVo_I = [dU(i) dV(i) dW(i)]';
    au = [OB_a.AX(i) OB_a.AY(i) OB_a.AZ(i)]';
    dVo_B(:,i) = R_BI*dVo_I;
    g_B(:,i) = R_BI*[0 0 9.812]';
    da(:,i) = R_BI*dVo_I - au - R_BI*[0 0 9.812]';
    p = OB_a.P(i); q = OB_a.Q(i); r = OB_a.R(i);
    p_dot = dP(i); q_dot = dQ(i); r_dot = dR(i);
    A = [q^2+r^2 -p*q+r_dot -p*r-q_dot;
         -p*q-r_dot p^2+r^2 -q*r+p_dot;
         -p*r+q_dot -q*r-p_dot p^2+q^2];
    datest(:,i) = A*rw';
%     rw(:,i) = A\da(:,i);
end

figure
subplot(3,1,1)
plot(OT_a.TIME,da(1,:)'); hold on;
plot(OT_a.TIME,datest(1,:)'); ylim([-10,10])
subplot(3,1,2)
plot(OT_a.TIME,da(2,:)'); hold on;
plot(OT_a.TIME,datest(2,:)');ylim([-10,10])
subplot(3,1,3)
plot(OT_a.TIME,da(3,:)'); hold on;
plot(OT_a.TIME,datest(3,:)'); ylim([-10,10])
% figure
% plot(OT_a.TIME,rw')
%%
i = 1;
n = 30;
x0 = lines(i).XData;
y0 = lines(i).YData;
x = compress_matrix(x0',n);
y = compress_matrix(y0',n);