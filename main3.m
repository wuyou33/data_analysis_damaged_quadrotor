%% For validating the model

clear all

for  index = 47; %b47 r54 81

addpath(genpath('E:\Data\_code_import_files'));
[OB_a,OT_a,WIND,PARA,take,DU] = import_data(index,1, 1, 1,1,1,'Q',0);

end

%%

%% preprocessing data

DU = [1, length(OT_a.TIME)];
% DU = [round(length(OT_a.TIME)/20), round(length(OT_a.TIME)*9.5/10)];
% DU = [10920 139600]; %RB n = 0;
% DU = [1 137900]; %LB n = 0
% DU = [24950 32450]; %#48 v = 1
% DU = [47680 51400]; %#48 v = 2
% DU = [62510 67440]; %#48 v = 3
% DU = [79030 82970]; %#48 v = 4
% DU = [91720 95970]; %#48 v = 5
% DU = [105700 109900]; %#48 v = 6
% DU = [120300 123900]; %#48 v = 7
% DU = [132600 139300]; %#48 v = 8

DU = [9729 13260]; %#47 v = 1
% DU = [28250 32120]; %#47 v = 2
% DU = [43380 46690]; %#47 v = 3
% DU = [57080 60620]; %#47 v = 4
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

% DU = [11710,12880]; %#135
du = DU(1):DU(2);

%% Use EKF to estimate rc


%%
addpath(genpath('E:\system identification\quadrotor identification\models'));
K = load('model_individual_type14');
N = length(OB_a.TIME);
l = PARA.l; b = PARA.b; m = PARA.mass/1000; Iv = PARA.Iv;
omega_mean = 1200^2;
omega1 = (OB_a.w1obs/60*2*pi).^2/omega_mean;
omega2 = (OB_a.w2obs/60*2*pi).^2/omega_mean;
omega3 = (OB_a.w3obs/60*2*pi).^2/omega_mean;
omega4 = (OB_a.w4obs/60*2*pi).^2/omega_mean;

% pqr = [OB_a.P,OB_a.Q,OB_a.R];
pqr = [butterworth(OB_a.p,4,15/256),butterworth(OB_a.q,4,15/256), butterworth(OB_a.r,4,15/256)];

V = [OT_a.U_air,OT_a.V_air,OT_a.W_air];

d = [ l,    -b, 0;
      l,     b, 0;
     -l,     b, 0;
     -l,    -b, 0];

u = V(:,1); 
v = V(:,2);
w = V(:,3);

p = pqr(:,1);
q = pqr(:,2);
r = pqr(:,3);

if mean(r) > 0
    rotation = 1;
else
    rotation = -1;
end

V1 = cross(pqr,repmat(d(1,:),size(V,1),1)) + V;
V2 = cross(pqr,repmat(d(2,:),size(V,1),1)) + V;
V3 = cross(pqr,repmat(d(3,:),size(V,1),1)) + V;
V4 = cross(pqr,repmat(d(4,:),size(V,1),1)) + V;

u1 = V1(:,1); u2 = V2(:,1); u3 = V3(:,1); u4 = V4(:,1);
v1 = V1(:,2); v2 = V2(:,2); v3 = V3(:,2); v4 = V4(:,2);
w1 = V1(:,3); w2 = V2(:,3); w3 = V3(:,3); w4 = V4(:,3);

SL = [1 -1 -1 1];
SM = [1 1 -1 -1];
SN = [-1 1 -1 1];

psi_mod = mod(OT_a.PSI+90,360); %% psi=0 point towards the nozzle
%% calculate force on the cg. (calculate r directly, assume only centrifugral force on IMU, failed)
acx = butterworth(OB_a.ax,4,5/256);
acy = butterworth(OB_a.ay,4,5/256);
acz = butterworth(OB_a.az,4,5/256);

gI = [0 0 9.8124]';
Ac = [acx,acy,acz];
rB = zeros(N,3);
for i = 1:N
G_centrifuge = [-(q(i).^2+r(i).^2), p(i).*q(i), p(i).*r(i);
                p(i).*q(i), -(r(i).*2+p(i).^2), q(i).*r(i);
                p(i).*r(i), q(i).*r(i), -(p(i).^2+q(i).^2)];
            
theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
        sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
        cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    

if (rank(G_centrifuge)==3)
rB(i,:) = (G_centrifuge\Ac(i,:)' + R_BI*gI)';   
else
    rB(i,:) = rB(i-1,:);
end
end

%% Calculate rc (estimate the displacement from optitrack centor to the cg by Least Square estimator)
% dVoB = fc + gB - \alpha x rc - \omega x (\omega x rc)

pqr_dot = [0 0 0; diff(pqr)*512];
dp = pqr_dot(:,1); dq = pqr_dot(:,2); dr = pqr_dot(:,3);
V_raw_inertia = [OT_a.vel_E(:,1),OT_a.vel_E(:,2),OT_a.vel_E(:,3)];
VoI = butterworth(V_raw_inertia,4,5.5/256);
dVoI = [0 0 0; diff(VoI)*512];
g = 9.8124;
fcx = butterworth(OB_a.ax,4,5/256)*g;
fcy = butterworth(OB_a.ay,4,5/256)*g;
fcz = butterworth(OB_a.az,4,5/256)*g;
% fc = [fcx, fcy, fcz];
fc = [fcx - fcx(1), fcy - fcy(1), fcz - (fcz(1)+9.8124)]; %  into account the bias

gB = zeros(size(VoI));
dVoB = zeros(size(VoI));
gI = [0 0 9.8124]';
rc = zeros(size(VoI));
rcxy = zeros(N,2);
GGr1 = zeros(size(VoI));
GGr2 = zeros(size(VoI));
GGr3 = zeros(size(VoI));
for i = DU(1):DU(2)
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
    
    R_IB = R_BI';   
    
    Gr1 = [0 -dr(i) dq(i); dr(i) 0 -dp(i); -dq(i) dp(i) 0];
    Gr2 = [-(q(i).^2+r(i).^2), p(i).*q(i), p(i).*r(i);
                p(i).*q(i), -(r(i).^2+p(i).^2), q(i).*r(i);
                p(i).*r(i), q(i).*r(i), -(p(i).^2+q(i).^2)];
    Gr = Gr1 + Gr2;
    
    GGr1(i,:) = [0 -dr(i) dq(i)] + [-(q(i).^2+r(i).^2), p(i).*q(i), p(i).*r(i)];
    GGr2(i,:) = [dr(i) 0 -dp(i)] + [p(i).*q(i), -(r(i).^2+p(i).^2), q(i).*r(i)];
    GGr3(i,:) = [-dq(i) dp(i) 0] + [p(i).*r(i), q(i).*r(i), -(p(i).^2+q(i).^2)];
    
%     GGr1(i,:) = [-(q(i).^2+r(i).^2), p(i).*q(i), p(i).*r(i)];
%     GGr2(i,:) = [p(i).*q(i), -(r(i).^2+p(i).^2), q(i).*r(i)];
%     GGr3(i,:) = [p(i).*r(i), q(i).*r(i), -(p(i).^2+q(i).^2)];
%     
    gB(i,:) = (R_BI*gI)';
    
    dVoB(i,:) = R_BI*dVoI(i,:)';
    rc(i,:) = Gr\(fc(i,:)'+gB(i,:)'-R_BI*dVoI(i,:)'); 
    
end

GGr = [GGr1(du,:); GGr2(du,:); GGr3(du,:)];
z = [fc(du,1)+gB(du,1)-dVoB(du,1);fc(du,2)+gB(du,2)-dVoB(du,2);fc(du,3)+gB(du,3)-dVoB(du,3)];
rcxy = lscov(GGr,z);

for i = DU(1):DU(2)
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
    
    R_IB = R_BI';   
    
    Gr1 = [0 -dr(i) dq(i); dr(i) 0 -dp(i); -dq(i) dp(i) 0];
    Gr2 = [-(q(i).^2+r(i).^2), p(i).*q(i), p(i).*r(i);
                p(i).*q(i), -(r(i).^2+p(i).^2), q(i).*r(i);
                p(i).*r(i), q(i).*r(i), -(p(i).^2+q(i).^2)];
    Gr = Gr1;
    
    Y(i,:) = (fc(i,:)' + gB(i,:)' - Gr*rcxy)';
end

figure
plot(dVoB(du,1),'b'); hold on;plot(dVoB(du,2) ,'r');plot(dVoB(du,3) ,'g');
 plot(Y(du,1),'b-.');plot(Y(du,2),'r-.');plot(Y(du,3),'g-.');
% YYY = fc(du,:)'+ gB(du,:)' - cross(pqr_dot',repmat(rcxy,1,N)) - cross(pqr,cross(pqr,repmat(rcxy,1,N)'));
%%
V_raw_inertia = [OT_a.vel_E(:,1),OT_a.vel_E(:,2),OT_a.vel_E(:,3)];
VoI = butterworth(V_raw_inertia,4,5.5/256);
VcI_filt = butterworth(VoI,4,2/256);
VoB = zeros(size(VoI));
VcB_filt = zeros(size(VcI_filt));

for i = 1:length(OT_a.TIME)
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
    
    R_IB = R_BI';
    
    VoB(i,:) = (R_BI * VoI(i,:)')';
    VcB_filt(i,:) = (R_BI * VcI_filt(i,:)')';
end

rpy = (VoB(:,1)-VcB_filt(:,1))./r;
rpx = (VcB_filt(:,2)-VoB(:,2))./r;
% rB = [rx, ry, zeros(size(rx))];

rpx = mean(rpx(6000:10000));
rpy = mean(rpy(6000:10000));
rpB = repmat([rpx, rpy, 0],length(OT_a.TIME),1);

pqr_filt = butterworth(pqr,4,15/256);
dVoI = [0,0,0;diff(VoI)]*512;
dVcI_filt = [0,0,0;diff(VcI_filt)]*512;

dVoB = zeros(size(dVoI));
dVcB2 = zeros(size(dVoI));
gB   = zeros(size(dVoI));

for i = 1:length(OT_a.TIME)
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
    
    R_IB = R_BI';
    
    dVoB(i,:) = (R_BI * dVoI(i,:)')';
    gB(i,:)   = (R_BI * [0 0 9.8124]')';
    
    % estimated dVc^B/dt by taking derivative of the filtered Vc^I; for comparison
    dVcB2(i,:)= (R_BI * dVcI_filt(i,:)')'; 
end

% (dVc)^B using method 1
dVcB = dVoB+ cross([0,0,0;diff(pqr_filt)]*512,rB) + cross(pqr_filt,cross(pqr_filt,rB));

% (dVc)^B/dt = d(Vc^B)/dt + omega x Vc^B; the outcome is the same as above.
VcB = VoB + cross(pqr,rB);
dVcB3 = [0,0,0;diff(VcB)]*512 + cross(pqr,VcB);

F = m*(dVcB - gB);
M = [0 0 0; diff(pqr)]*Iv;  %% warning!! the Iv is for original bb2 instead of the lightweighted one.

%%
%% WARNING !! THE ROTATION CENTER IS NOT CG! IMAGINE THE WOBBLING MOTION
%% THE VcB is actually the velocity of the spinning center instead of the CG.

%% model verification

Fx = F(:,1);
Fy = F(:,2);
Fz = F(:,3);
Mx = M(:,1);
My = M(:,2);
Mz = M(:,3);

load('E:\system identification\quadrotor identification\models\model_individual_type14');

if rotation == 1
    fx = kx*u1.*omega1 + kx*u2.*omega2 + kx*u3.*omega3; % + kx*u4.*omega4;
    fy = ky*v1.*omega1 + ky*v2.*omega2 + ky*v3.*omega3; % + ky*v4.*omega4;   
else
    fx = kx*u1.*omega1 + kx*u2.*omega2 + kx*u4.*omega4;
    fy = ky*v1.*omega1 + ky*v2.*omega2 + ky*v4.*omega4;    
end


return;
%% PLOT
plot_psi_vs_f = 1;
if plot_psi_vs_f
   figure
   subplot(2,1,1)
   plot(psi_mod(du),Fx(du),'.'); hold on;
   plot(psi_mod(du),fx(du),'.'); ylabel('fx'); legend('mea','model');
   subplot(2,1,2)
   plot(psi_mod(du),Fy(du),'.'); hold on;
   plot(psi_mod(du),fy(du),'.'); ylabel('fy'); legend('mea','model');
end