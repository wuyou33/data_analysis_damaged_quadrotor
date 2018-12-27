
clear all

addpath('E:\surfdrive\DATA\func_common');
addpath('kf_func');
addpath(genpath('E:\surfdrive\DATA\_data_cyberzoo'));

load('take.mat');
take.name = ['2017-10-19 10.38.15 AM'];
take.rigid_body = 'bebop2_sihao';

filt_OT = 1;
OT_orig         = import_optitrack_data(take.name,take.name_calib,{take.rigid_body},filt_OT);
OB              = import_onboard_data(take.name);
% OB              = import_onboard_data2(take.name);
OT_filt         = OT_orig;

% OT_a = OT_orig;
% % OB_a = OB;
%% import wind speed
T_wind = table2array(readtable(['windspeed ' take.name '.txt'],'Delimiter','\t'));
time_wind = T_wind(:,1);
speed_wind = T_wind(:,5);
time_wind(time_wind==0) = []; speed_wind(speed_wind==0) = [];

take_date = str2double(take.name(9:10));
take_hour = str2double(take.name(12:13));
take_minute = str2double(take.name(15:16));
take_sec = str2double(take.name(18:19));
take_AMPM = take.name(21:22);
if take_AMPM == 'PM'
    take_hour = take_hour + 12;
end
% relevate time to 2017-11-08 13.00.00 PM
take_time = (take_date - 8)*3600*24 + (take_hour - 13)*3600 + take_minute*60 + take_sec; 


%% resample and finddelay
OT = struct;
allnames = fieldnames(OT_filt);
OB_discard = find(OB.TIME > OT_filt.TIME(end));
time_OB = OB.TIME;
time_OB(OB_discard) = [];
for i = 1:length(fieldnames(OT_filt))
    ct_orig = OT_filt.(allnames{i});
    ct_orig = interpNan(ct_orig);
    OT.(allnames{i}) = interp1(OT_filt.TIME,ct_orig,time_OB,'spline');
end
OT.TIME = time_OB;

[OT_a,OB_a,Delay] = align_signals(OT,OB,'phi');

time_wind_a = time_wind - (take_time + Delay/512);
speed_wind_a = speed_wind;
speed_wind_a(time_wind_a<0) = [];
time_wind_a(time_wind_a<0) = [];
%% Add wind speed to OT_a
speed_wind_resemple = interp1(time_wind_a,speed_wind_a,OT_a.TIME,'spline');
VY_wind = OT_a.VY - speed_wind_resemple;
%% preprocessing data
wsum = OB_a.w1obs+OB_a.w2obs+OB_a.w3obs;

Vz_filt = butterworth(OT_a.VZ,4,4/256);
azi_filt = derivative(Vz_filt,OB_a.TIME);
psi_mod = mod(OT_a.PSI,360)-180;

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

figure
subplot(2,1,1)
plot(OB_a.TIME,OB_a.h1); hold on; grid on; title(['h1,h2',' ',take.name]);
% xlim([22,24]);
subplot(2,1,2)
plot(OB_a.TIME,OB_a.h2); grid on;
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

% figure
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
plot(OB_a.TIME,OB_a.w4ref); hold on; grid on; title(['rotor speed',' ',take.name]);
plot(OB_a.TIME,OB_a.w4obs);
% xlim([22,24]);

figure
subplot(3,1,1)
plot(OB_a.TIME,OB_a.p_des); hold on; grid on; title(['pq_{des} and pq',' ',take.name]);
plot(OB_a.TIME,OB_a.p);
% xlim([22,24]);
% plot(OT_a.TIME,OT_a.P);
subplot(3,1,2)
plot(OB_a.TIME,OB_a.q_des); hold on; grid on;
plot(OB_a.TIME,OB_a.q);
% xlim([22,24]);
% plot(OT_a.TIME,OT_a.Q);
subplot(3,1,3)
plot(OB_a.TIME,OB_a.r,'r'); hold on; grid on;

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
%% 
DU = [1, length(OT_a.TIME)];
DU = [round(length(OT_a.TIME)/20), round(length(OT_a.TIME)*9.5/10)];
% DU = [10920 139600]; %RB n = 0;
% DU = [1 137900]; %LB n = 0
fc = 2;
du = DU(1):DU(2);

% calculate acceleartion from external sensor (Optitrack)
VX_filt = butterworth(OT_a.VX,4,fc/256);
VY_filt = butterworth(OT_a.VY,4,fc/256);
VZ_filt = butterworth(OT_a.VZ,4,fc/256);
VX_filt = butterworth(VX_filt,4,fc/256);
VY_filt = butterworth(VY_filt,4,fc/256);
VZ_filt = butterworth(VZ_filt,4,fc/256);
VY_wind_filt = butterworth(VY_wind,4,fc/256);

dVX = derivative(VX_filt,OT_a.TIME);
dVY = derivative(VY_filt,OT_a.TIME);
dVZ = derivative(VZ_filt,OT_a.TIME);

% project n on Inertia frame; project velocity and dV on Body frame. 
na_I = zeros(length(OT_a.TIME),2); 
% n = [0.2 0.2 -1.0]';
n   = [0.0 0.0 -1.0]'; %primary axis
n_I = zeros(length(OT_a.TIME),3);
Ab  = zeros(length(OT_a.TIME),3); 
Vb  = zeros(length(OT_a.TIME),3); 
for i = 1:length(OT_a.TIME)
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
    
    R_IB = [cos(psi)*cos(theta) , cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);
        sin(psi)*cos(theta) , sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);
        -sin(theta)          , cos(theta)*sin(phi)                          , cos(theta)*cos(phi)                          ];
    
    Ab(i,:) = [dVX(i) dVY(i) dVZ(i)]*R_BI';
    Vb(i,:) = [VX_filt(i) VY_filt(i) VZ_filt(i)]*R_BI';
    n_I(i,:) =  n'*R_IB';
end

% calculate Stability frame. Project Velocity and dV on Stability frame.
na_I(:,1) = butterworth(n_I(:,1),4,fc/512);
na_I(:,2) = butterworth(n_I(:,2),4,fc/512);

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
    Vs(i,:) = [VX_filt(i) VY_wind_filt(i) VZ_filt(i)] * R_SI';
end

% plot figures

figure
plot(OT_a.TIME,phi_s);
hold on
plot(OT_a.TIME,theta_s);
xlabel('TIME'); legend('phi_s','theta_s');

figure
% subplot(2,2,1);
plot(Vs(du,2),As(du,1)); xlabel('Vy_{air}'); ylabel('As');
hold on;
plot(Vs(du,2),As(du,2));
plot(Vs(du,2),As(du,3));
title(take.name);
% figure
% % subplot(2,2,1);
% plot(Vs(du,1),As(du,1)./wsum(du));
% hold on;
% plot(Vs(du,2),As(du,2)./wsum(du));

% figure
% subplot(2,1,1)
% plot(Vb(du,1),Ab(du,1)); grid on; hold on; ylabel('acc^b_x'); xlabel('u^b');
% subplot(2,1,2)
% plot(Vb(du,2),Ab(du,2)); grid on; hold on; ylabel('acc^b_y'); xlabel('v^b');

return;
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