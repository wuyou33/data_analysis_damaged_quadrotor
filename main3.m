%% 
% This script estimate the forces and moments on the quadrotor by using the
% measurements from the onboard sensor. Centrifugral forces / displacement
% from c.g to the IMU are taken into account. The measuremens are stored in
% the ./data folder. The scripts in the ./system identification folder will
% use these data to identify forces and moment models.
%%
import = 1;

%%
if import 
    
    clear all;
    
    for  index = 9
        
    addpath(genpath('E:\surfdrive\DATA\_code_import_files'));
    file_index = 'E:\surfdrive\DATA\_data_OJF_june_2018\SRFDRF\Log_OJF_June';
%     file_index = 'E:\surfdrive\DATA\_data_OJF_may_2018\Log_OJF_May';
%     file_index = 'E:\surfdrive\DATA\_data_OJF_november_2017\Log_OJF_Nov';
    [OB_a,OT_a,WIND,PARA,take,DU] = import_data(index,file_index,0,'phi_ot', 1, 0, 1, 0);
    end
    
    import = 1;
end


validate_force = 1; 
%%
DU(1)
%% select data that

% DU = [1, length(OT_a.TIME)];
% DU = [round(length(OT_a.TIME)/30), round(length(OT_a.TIME)*9.0/10)];

% DRF 2#4 removed, r>0
DU = [25120 33700]; %#9 V = 2
% DU = [46430 54350]; %#9 V = 4
% DU = [65850 72360]; %#9 V = 6
% DU = [83100 91680]; %#9 V = 8
% DU = [9877 91680];

% DU = [23480 30410]; %#17 V = 2
% DU = [46080 50270]; %#17 V = 4
% DU = [65220 70390]; %#17 V = 6
% DU = [80610 86200]; %#17 V = 8
% DU = [8445,86200];

% DU = [27030 32830]; %#18 V = 2
% DU = [44940 51590]; %#18 V = 4
% DU = [62130 69970]; %#18 V = 6
% DU = [81820 88380]; %#18 V = 8
% DU = [6865 88380];

% DU = [35590 43100]; %#75 V = 2 DRF May
% DU = [60840 74030]; %#75 V = 4
% DU = [94620 96290]; %#75 V = 6
% DU = [114700 119800]; %#75 V = 8
% DU = [35590 119800];

% DRF 1#3 removed, r<0
% DU = [36850 43360]; %#7 V = 2
% DU = [58860 63110]; %#7 V = 4
% DU = [73680 79910]; %#7 V = 6
% DU = [90800 96460]; %#7 V = 8
% DU = [23650 96460];

% DU = [27950 34420]; %#19 V = 2
% DU = [43350 51370]; %#19 V = 4
% DU = [63870 69490]; %#19 V = 6
% DU = [80510 88220]; %#19 V = 8
% DU = [12350 88220];

% DU = [23850 28950]; %#20 V = 2
% DU = [41890 49450]; %#20 V = 4
% DU = [68020 71420]; %#20 V = 6
% DU = [79690 87430]; %#20 V = 8
% DU = [7174 87430];


% DU = [29110 37430]; %#32 V = 2 SRF r<0 LF removed
% DU = [50000 57940]; %#32 V = 4
% DU = [68230 73580]; %#32 V = 6
% DU = [83230 84740]; %#32 V = 8
% DU = [29110 84740];

% DU = [22600 30440]; %#35 V = 2 SRF r>0 rF removed
% DU = [39420 50420]; %#35 V = 4
% DU = [52710 53420]; %#35 V = 6
% DU = [22600 53420];

% DU = [30700 42220]; %#23 v = 2 SRF r>0 lb removed
% DU = [57670 65810]; %#23 v = 4
% DU = [77980 83230]; %#23 v = 6
% DU = [98090 100200]; %#23 v = 8
% DU = [3891 100200];

% DU = [31500 39470]; %#25 v = 2 SRF r<0 rb removed
% DU = [53140 58640]; %#25 v = 4
% DU = [67680 75360]; %#25 v = 6
% DU = [85470 92750]; %#25 v = 8
% DU = [31500 92750];


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
% DU = [57080 60620]; %#47 v = 4
% DU = [73450 76260]; %#47 v = 5
% DU = [83860 87020]; %#47 v = 6
% DU = [93300 96690]; %#47 v = 7
DU = [105400 106000]; %#47 v = 8

% DU = [24380 90710]; %#92 2-10 Norminal

% DU = [1 63210]; %#93 4-10 Norminal beta = -90;

% DU = [1 161100]; %#122 up_down nominal;

du = DU(1):DU(2);

%% Use EKF to estimate rc


%%
addpath(genpath('E:\system identification\quadrotor identification\models'));
K = load('model_individual');
N = length(OB_a.TIME);
l = PARA.l; b = PARA.b; m = PARA.mass/1000; Iv = PARA.Iv; R = PARA.R;
Area = pi*R^2;
rho = 1.225;

omega_mean = 1200^2;

omega1_b = butterworth(OB_a.w1obs/60*2*pi,4,15/256);
omega2_b = butterworth(OB_a.w2obs/60*2*pi,4,15/256);
omega3_b = butterworth(OB_a.w3obs/60*2*pi,4,15/256);
omega4_b = butterworth(OB_a.w4obs/60*2*pi,4,15/256);

omega1 = omega1_b - butterworth(OB_a.r,4,15/256);
omega2 = omega2_b + butterworth(OB_a.r,4,15/256);
omega3 = omega3_b - butterworth(OB_a.r,4,15/256);
omega4 = omega4_b + butterworth(OB_a.r,4,15/256);

if mean(OB_a.w1ref)<2800
    omega1_b = 0*omega1_b;
    omega1 = 0*omega1_b;
end
if mean(OB_a.w2ref)<2800
    omega2_b = 0*omega2_b;
    omega2 = 0*omega2_b;
end
if mean(OB_a.w3ref)<2800
    omega3_b = 0*omega3_b;
    omega3 = 0*omega3_b;
end
if mean(OB_a.w4ref)<2800
    omega4_b = 0*omega4_b;
    omega4 = 0*omega4_b;
end

omega1_dot = [0;diff(omega1_b)]*512;
omega2_dot = [0;diff(omega2_b)]*512;
omega3_dot = [0;diff(omega3_b)]*512;
omega4_dot = [0;diff(omega4_b)]*512;

omega1_bar = omega1.^2/omega_mean;
omega2_bar = omega2.^2/omega_mean;
omega3_bar = omega3.^2/omega_mean;
omega4_bar = omega4.^2/omega_mean;


% pqr = [OB_a.P,OB_a.Q,OB_a.R];
pqr = [butterworth(OB_a.p,4,15/256),butterworth(OB_a.q,4,15/256), butterworth(OB_a.r,4,15/256)];
% pqr = [butterworth(OT_a.P,4,15/256),butterworth(OT_a.Q,4,15/256), butterworth(OT_a.R,4,15/256)];

V = [OT_a.U_air,OT_a.V_air,OT_a.W_air];

u = V(:,1); 
v = V(:,2);
w = V(:,3);

va = sqrt(u.^2+v.^2+w.^2);

p = pqr(:,1);
q = pqr(:,2);
r = pqr(:,3);

if mean(r) > 0
    rotation = 1; %#2 or #4 removed
else
    rotation = -1; %#1 or #3 removed
end

switch take.configuration
    case 'drf'
           Iv(1,1) = Iv(1,1) + b^2*0.005*2;
           Iv(2,2) = Iv(2,2) + l^2*0.005*2;
           Iv(3,3) = Iv(3,3) + (b^2+l^2)*0.005*2;
    case 'srf'
           Iv(1,1) = Iv(1,1) + b^2*0.005*3;
           Iv(2,2) = Iv(2,2) + l^2*0.005*3;
           Iv(3,3) = Iv(3,3) + (b^2+l^2)*0.005*3;
    otherwise
end
Ip = PARA.Ip(3,3);

SL = [1 -1 -1 1];
SM = [1 1 -1 -1];
SN = [-1 1 -1 1];

psi_mod = mod(OT_a.PSI+90,360); %% psi=0 point towards the nozzle
psi_mod2 = mod(OT_a.PSI+90,720);

%% Calculate rc (estimate the displacement from optitrack centor to the IMU by Least Square estimator) (IN THSI CASE, ASSUME IMU is located at c.g)
% dVoB = fc + gB - \alpha x rc - \omega x (\omega x rc)

pqr_dot = [0 0 0; diff(pqr)*512];
dp = pqr_dot(:,1); dq = pqr_dot(:,2); dr = pqr_dot(:,3);
V_raw_inertia = [OT_a.vel_E(:,1),OT_a.vel_E(:,2),OT_a.vel_E(:,3)];
VoI = butterworth(V_raw_inertia,4,15/256);
% VoI = V_raw_inertia;
dVoI = [0 0 0; diff(VoI)*512];
g = 9.8124;

if strcmp(take.drone,'bebop2_guido')
    dax=-0.198937;
    day=-0.118300;
    daz=-0.314115;
else
    dax = 0; day = 0; daz = 0;
end

if strcmp(file_index, 'E:\surfdrive\DATA\_data_OJF_november_2017\Log_OJF_Nov')
    fcx = butterworth(OB_a.ax,4,15/256) - dax;
    fcy = butterworth(OB_a.ay,4,15/256) - day;
    fcz = butterworth(OB_a.az,4,15/256) - daz;
elseif strcmp(file_index, 'E:\surfdrive\DATA\_data_OJF_may_2018\Log_OJF_May') && strcmp(take.configuration,'drf')
    fcx = butterworth(OB_a.ax,4,15/256) - dax;
    fcy = butterworth(OB_a.ay,4,15/256) - day;
    fcz = butterworth(OB_a.az,4,15/256) - daz;    
else
    fcx = butterworth(OB_a.ax,4,15/256)*g - dax;
    fcy = butterworth(OB_a.ay,4,15/256)*g - day;
    fcz = butterworth(OB_a.az,4,15/256)*g - daz;
end
% fcx = butterworth(OB_a.ax,4,15/256);
% fcy = butterworth(OB_a.ay,4,15/256);
% fcz = butterworth(OB_a.az,4,15/256);    

fc = [fcx, fcy, fcz];
% fc = [fcx - fcx(1), fcy - fcy(1), fcz - (fcz(1)+9.8124)]; %  into account the bias

gB = zeros(size(VoI));
dVoB = zeros(size(VoI));
gI = [0 0 9.8124]';
rcxy = zeros(N,2);
GGr1 = zeros(size(VoI));
GGr2 = zeros(size(VoI));
GGr3 = zeros(size(VoI));
Gr = zeros(3,3,length(VoI));
for i = 1:N
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
    
    R_IB = R_BI';   
    
    Gr1 = [0 -dr(i) dq(i); dr(i) 0 -dp(i); -dq(i) dp(i) 0];
    Gr2 = [-(q(i).^2+r(i).^2), p(i).*q(i), p(i).*r(i);
                p(i).*q(i), -(r(i).^2+p(i).^2), q(i).*r(i);
                p(i).*r(i), q(i).*r(i), -(p(i).^2+q(i).^2)];
    Gr(:,:,i) = Gr1 + Gr2;
    
    GGr1(i,:) = [0 -dr(i) dq(i)] + [-(q(i).^2+r(i).^2), p(i).*q(i), p(i).*r(i)];
    GGr2(i,:) = [dr(i) 0 -dp(i)] + [p(i).*q(i), -(r(i).^2+p(i).^2), q(i).*r(i)];
    GGr3(i,:) = [-dq(i) dp(i) 0] + [p(i).*r(i), q(i).*r(i), -(p(i).^2+q(i).^2)];
    
%     GGr1(i,:) = [-(q(i).^2+r(i).^2), p(i).*q(i), p(i).*r(i)];
%     GGr2(i,:) = [p(i).*q(i), -(r(i).^2+p(i).^2), q(i).*r(i)];
%     GGr3(i,:) = [p(i).*r(i), q(i).*r(i), -(p(i).^2+q(i).^2)];
%     
    
    gB(i,:) = (R_BI*gI)';
    
    dVoB(i,:) = R_BI*dVoI(i,:)';
%     rc(i,:) = Gr\(fc(i,:)'+gB(i,:)'-R_BI*dVoI(i,:)'); 
    
end

GGr = [GGr1(du,:); GGr2(du,:); GGr3(du,:)];
GGrflat = [GGr1(du,:); GGr2(du,:)];
zrc = [fc(du,1)+gB(du,1)-dVoB(du,1);fc(du,2)+gB(du,2)-dVoB(du,2);fc(du,3)+gB(du,3)-dVoB(du,3)];
zra_flat = [fc(du,1); fc(du,2)];
% zra = [fc(du,1); fc(du,2)]
if import
    rc = lscov(GGr,zrc);
end

if import  % distance from c.g to imu
%     ra = lscov(GGr,zra);
    ra = lscov(GGrflat,zra_flat);
end

% validate the correctness of rc
if import == 1
VcB = zeros(size(VoI));
VoB = zeros(size(VoI));

% for i = DU(1):DU(2)
for i = 1:length(VoI)
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
    
    R_IB = R_BI';      
    
    Y(i,:) = (fc(i,:)' + gB(i,:)' - Gr(:,:,i)*rc)';
    VoB(i,:) = (R_BI*VoI(i,:)')';
    VcB(i,:) = VoB(i,:) + cross(pqr(i,:),rc);

end

figure
plot(dVoB(du,1),'b'); hold on;plot(dVoB(du,2) ,'r');plot(dVoB(du,3) ,'g');
plot(Y(du,1),'b-.');plot(Y(du,2),'r-.');plot(Y(du,3),'g-.');
% YYY = fc(du,:)'+ gB(du,:)' - cross(pqr_dot',repmat(rcxy,1,N)) - cross(pqr,cross(pqr,repmat(rcxy,1,N)'));

end
%% Plot resultant force on the c.g. with the information of accelerometer and ra.

f = zeros(size(VoI));
for i = 1:N
    
    f(i,:) = fc(i,:) + (-Gr(:,:,i)*ra)';

end

figure
subplot(3,1,1)
plot(fc(du,1)); hold on;
plot(f(du,1)); ylabel('a_x [m/s^2]')
subplot(3,1,2)
plot(fc(du,2)); hold on;
plot(f(du,2)); ylabel('a_y [m/s^2]')
subplot(3,1,3)
plot(fc(du,3)); hold on;
plot(f(du,3)); ylabel('a_z [m/s^2]')
legend('mea','corrected')

if validate_force == 1
F = m*f;
M = [0 0 0; diff(pqr)*512]*Iv + cross(pqr',Iv*pqr')' + ...
    [Ip*q.*(-omega1_b+omega2_b-omega3_b+omega4_b), -Ip*p.*(-omega1_b+omega2_b-omega3_b+omega4_b),zeros(size(p))];
M(:,3) = M(:,3) + Ip * (-omega1_dot+omega2_dot-omega3_dot+omega4_dot);
% M = [0 0 0; diff(pqr)*512]*Iv;
% Fx = F(:,1);
% Fy = F(:,2);
% Fz = F(:,3);
Fx = f(:,1)*m;
Fy = f(:,2)*m;
Fz = f(:,3)*m;
Mx = M(:,1);
My = M(:,2);
Mz = M(:,3);


%%
d = [ l + rc(1),    -b + rc(2), 0 + rc(3);
      l + rc(1),     b + rc(2), 0 + rc(3);
     -l + rc(1),     b + rc(2), 0 + rc(3);
     -l + rc(1),    -b + rc(2), 0 + rc(3)];
 
V1 = cross(pqr,repmat(d(1,:),size(V,1),1)) + V;
V2 = cross(pqr,repmat(d(2,:),size(V,1),1)) + V;
V3 = cross(pqr,repmat(d(3,:),size(V,1),1)) + V;
V4 = cross(pqr,repmat(d(4,:),size(V,1),1)) + V;

u1 = V1(:,1); u2 = V2(:,1); u3 = V3(:,1); u4 = V4(:,1);
v1 = V1(:,2); v2 = V2(:,2); v3 = V3(:,2); v4 = V4(:,2);
w1 = V1(:,3); w2 = V2(:,3); w3 = V3(:,3); w4 = V4(:,3);
va1 = sqrt(u1.^2+v1.^2+w1.^2);
va2 = sqrt(u2.^2+v2.^2+w2.^2);
va3 = sqrt(u3.^2+v3.^2+w3.^2);
va4 = sqrt(u4.^2+v4.^2+w4.^2);

load('E:\system identification\quadrotor identification\models\Cz_airframe_BB2');
% load('E:\system identification\Rotor_Data\rotor_static_drag_model_BB2.mat');
load('E:\system identification\damaged_model_identification\Rotor_Data\rotor_torque_model_BB2.mat');
load('E:\system identification\damaged_model_identification\Ct_Cq_model\Ct_model_BB2.mat');
load('E:\system identification\damaged_model_identification\Ct_Cq_model\Cq_model_BB2.mat');

alpha = atan(w./(sqrt(u.^2+v.^2)))*57.3;
alpha1 = atan(w1./sqrt(u1.^2+v1.^2))*57.3;
alpha2 = atan(w2./sqrt(u2.^2+v2.^2))*57.3;
alpha3 = atan(w3./sqrt(u3.^2+v3.^2))*57.3;
alpha4 = atan(w4./sqrt(u4.^2+v4.^2))*57.3;

fz_airframe = (u.^2 + v.^2 + w.^2).* interp1(AoA_airframe, Cz_airframe, alpha ,'spline');
% fz_rotor_static = interp2(Vel_dense_grid,AoA_dense_grid,thrust_dense_vali_rotor_net,...
%                             va1,alpha1,'spline') + ...
%                   interp2(Vel_dense_grid,AoA_dense_grid,thrust_dense_vali_rotor_net,...
%                             va2,alpha2,'spline') + ...
%                   interp2(Vel_dense_grid,AoA_dense_grid,thrust_dense_vali_rotor_net,...
%                             va3,alpha3,'spline') + ...
%                   interp2(Vel_dense_grid,AoA_dense_grid,thrust_dense_vali_rotor_net,...
%                             va4,alpha4,'spline');

% fz_airframe = 0;

fx = K.kx*u1.*omega1_bar + K.kx*u2.*omega2_bar + K.kx*u3.*omega3_bar + K.kx*u4.*omega4_bar;
fy = K.ky*v1.*omega1_bar + K.ky*v2.*omega2_bar + K.ky*v3.*omega3_bar + K.ky*v4.*omega4_bar;


% fz from static rotor data.
vv1 = va1./omega1/R;
vv2 = va2./omega2/R;
vv3 = va3./omega3/R;
vv4 = va4./omega4/R;

i_inf1 = intersect(find(vv1 ~= inf),find(vv1 ~= -inf));
i_inf2 = intersect(find(vv2 ~= inf),find(vv2 ~= -inf));
i_inf3 = intersect(find(vv3 ~= inf),find(vv3 ~= -inf));
i_inf4 = intersect(find(vv4 ~= inf),find(vv4 ~= -inf));

Ct1 = zeros(size(vv1)); Ct2 = zeros(size(vv2)); Ct3 = zeros(size(vv3)); Ct4 = zeros(size(vv4));
Ct1(i_inf1) = Ct_model_BB2_surf(alpha1(i_inf1),vv1(i_inf1));
Ct2(i_inf2) = Ct_model_BB2_surf(alpha2(i_inf2),vv2(i_inf2));
Ct3(i_inf3) = Ct_model_BB2_surf(alpha3(i_inf3),vv3(i_inf3));
Ct4(i_inf4) = Ct_model_BB2_surf(alpha4(i_inf4),vv4(i_inf4));

T1 = Ct1.*(omega1*R).^2*Area*rho;
T2 = Ct2.*(omega2*R).^2*Area*rho;
T3 = Ct3.*(omega3*R).^2*Area*rho;
T4 = Ct4.*(omega4*R).^2*Area*rho;

fz = -(T1 + T2 + T3+ T4) + fz_airframe;

% fz = -(P33(u1,v1,w1,omega1_bar)*K.kt + P33(u2,v2,w2,omega2_bar)*K.kt + P33(u3,v3,w3,omega3_bar)*K.kt + P33(u4,v4,w4,omega4_bar)*K.kt) ...
%      + fz_airframe + fz_rotor_static;
%  
% fz = -((P33(u1,v1,w1)+P33(u2,v2,w2)+P33(u3,v3,w3)+P33(u4,v4,w4))*K.kt0 + P33(u1,v1,w1,omega1_bar)*K.kt + P33(u2,v2,w2,omega2_bar)*K.kt + P33(u3,v3,w3,omega3_bar)*K.kt + P33(u4,v4,w4,omega4_bar)*K.kt);

% take into account dT effect
% dCt_model = load('E:\system identification\thrust_model\dCt_model_BB2.mat');
mu1 = sqrt(u1.^2+v1.^2)./(omega1.*R); mu1(abs(mu1)==inf)=0;
mu2 = sqrt(u2.^2+v2.^2)./(omega2.*R); mu2(abs(mu2)==inf)=0;
mu3 = sqrt(u3.^2+v3.^2)./(omega3.*R); mu3(abs(mu3)==inf)=0;
mu4 = sqrt(u4.^2+v4.^2)./(omega4.*R); mu4(abs(mu4)==inf)=0;

lc1 = w1./(omega1.*R); lc1(abs(lc1)==inf)=0;
lc2 = w2./(omega2.*R); lc2(abs(lc2)==inf)=0; 
lc3 = w3./(omega3.*R); lc3(abs(lc3)==inf)=0;
lc4 = w4./(omega4.*R); lc4(abs(lc4)==inf)=0;

psi_h1 = psi_mod - 417; psi_h1 = psi_h1/57.3;
psi_h2 = psi_mod - 308; psi_h2 = psi_h2/57.3;
psi_h3 = psi_mod - 232; psi_h3 = psi_h3/57.3;
psi_h4 = psi_mod - 128; psi_h4 = psi_h4/57.3;

psi_h1 = mod(psi_h1,2*pi);
psi_h2 = mod(psi_h2,2*pi);
psi_h3 = mod(psi_h3,2*pi);
psi_h4 = mod(psi_h4,2*pi);

dynhead1 = rho.*omega1.^2.*R^2;
dynhead2 = rho.*omega2.^2.*R^2;
dynhead3 = rho.*omega3.^2.*R^2;
dynhead4 = rho.*omega4.^2.*R^2;

dT1 = 0; dT2 = 0; dT3 = 0; dT4 = 0;
%%
mx = SL(1)*b*(T1+dT1) + SL(2)*b*(T2+dT2) + SL(3)*b*(T3+dT3) + SL(4)*b*(T4+dT4);
my = SM(1)*l*(T1+dT1) + SM(2)*l*(T2+dT2) + SM(3)*l*(T3+dT3) + SM(4)*l*(T4+dT4);

% mx = (SL(1)*b*P33(u1,v1,w1,omega1_bar)+SL(2)*b*P33(u2,v2,w2,omega2_bar)+SL(3)*b*P33(u3,v3,w3,omega3_bar))*K.kt ...
%      + (P33(SN(1)*u1,v1,w1,omega1_bar)+P33(SN(2)*u2,v2,w2,omega2_bar)+P33(SN(3)*u3,v3,w3,omega3_bar))*K.kl;
% my = SM(1)*l*P33(u1,v1,w1,omega1_bar)*K.kt +SM(2)*l*P33(u2,v2,w2,omega2_bar)*K.kt+SM(3)*l*P33(u3,v3,w3,omega3_bar)*K.kt ...
%      + (P33(u1,SN(1)*v1,w1,omega1_bar)+P33(u2,SN(2)*v2,w2,omega2_bar)+P33(u3,SN(3)*v3,w3,omega3_bar))*K.km;
% mz = b*(u1.*omega1_bar-u2.*omega2_bar-u3.*omega3_bar+u4.*omega4_bar)*K.kx + K.ky*l*(v1.*omega1_bar+v2.*omega2_bar-v3.*omega3_bar-v4.*omega4_bar)...
%      + (SN(1)*P33(u1,v1,w1,omega1_bar)+SN(2)*P33(u2,v2,w2,omega2_bar)+SN(3)*P33(u3,v3,w3,omega3_bar)+SN(4)*P33(u4,v4,w4,omega4_bar))*K.kn ...
%      + P33(u,v,w) * K.kn0;

% mz from the single rotor data.
% Cq01 = 0; Cq02 = 0; Cq03=0; Cq04=0;
% if mean(omega1_bar ~= 0)
% Cq01 = interp2(Vel_sparse_grid,AoA_sparse_grid,Cq0_rotor_torque_model_sparse,...
%                             va1,alpha1,'spline');
% end
% if mean(omega2_bar ~= 0)
% Cq02 = interp2(Vel_sparse_grid,AoA_sparse_grid,Cq0_rotor_torque_model_sparse,...
%                             va2,alpha2,'spline');
% end
% if mean(omega3_bar ~= 0)
% Cq03 = interp2(Vel_sparse_grid,AoA_sparse_grid,Cq0_rotor_torque_model_sparse,...
%                             va3,alpha3,'spline');
% end
% if mean(omega4_bar ~= 0)
% Cq04 = interp2(Vel_sparse_grid,AoA_sparse_grid,Cq0_rotor_torque_model_sparse,...
%                             va4,alpha4,'spline');
% end                        
% Cq1 = interp2(Vel_sparse_grid,AoA_sparse_grid,Cq_rotor_torque_model_sparse,...
%                             va1,alpha1,'spline');
% Cq2 = interp2(Vel_sparse_grid,AoA_sparse_grid,Cq_rotor_torque_model_sparse,...
%                             va2,alpha2,'spline');
% Cq3 = interp2(Vel_sparse_grid,AoA_sparse_grid,Cq_rotor_torque_model_sparse,...
%                             va3,alpha3,'spline');
% Cq4 = interp2(Vel_sparse_grid,AoA_sparse_grid,Cq_rotor_torque_model_sparse,...
%                             va4,alpha4,'spline');
% mz = SN(1)*(Cq01+Cq1.*omega1_bar*omega_mean) + SN(2)*(Cq02+Cq2.*omega2_bar*omega_mean) + ...
%      SN(3)*(Cq03+Cq3.*omega3_bar*omega_mean) + SN(4)*(Cq04+Cq4.*omega4_bar*omega_mean);

% mz from the dimensionless single rotor data
Cq1 = zeros(size(vv1)); Cq2 = zeros(size(vv2)); Cq3 = zeros(size(vv3)); Cq4 = zeros(size(vv4));
Cq1(i_inf1) = Cq_model_BB2_surf(alpha1(i_inf1),vv1(i_inf1));
Cq2(i_inf2) = Cq_model_BB2_surf(alpha2(i_inf2),vv2(i_inf2));
Cq3(i_inf3) = Cq_model_BB2_surf(alpha3(i_inf3),vv3(i_inf3));
Cq4(i_inf4) = Cq_model_BB2_surf(alpha4(i_inf4),vv4(i_inf4));

Q1 = -Cq1.*(omega1*R).^2*Area*rho;
Q2 =  Cq2.*(omega2*R).^2*Area*rho;
Q3 = -Cq3.*(omega3*R).^2*Area*rho;
Q4 =  Cq4.*(omega4*R).^2*Area*rho;

mz = Q1 + Q2 + Q3 + Q4;
 %%

% figure
% subplot(3,1,1)
% plot(psi_mod(du),Fx(du),'b.'); hold on;
% plot(psi_mod(du),fx(du),'r.'); ylabel('fx'); legend('mea','model');
% subplot(3,1,2)
% plot(psi_mod(du),Fy(du),'b.'); hold on;
% plot(psi_mod(du),fy(du),'r.'); ylabel('fy'); legend('mea','model');
% subplot(3,1,3)
% plot(psi_mod(du),Fz(du),'b.'); hold on; ylabel('fz');
% plot(psi_mod(du),fz(du),'r.'); hold on; legend('mea','model');


figure
subplot(2,2,1)
plot(psi_mod(du),Mx(du),'.'); hold on;
plot(psi_mod(du),mx(du),'.');
subplot(2,2,2)
plot(Mx(du)); hold on;
plot(mx(du));
subplot(2,2,3)
plot(psi_mod(du),My(du),'.'); hold on;
plot(psi_mod(du),my(du),'.');
subplot(2,2,4)
plot(My(du)); hold on;
plot(my(du));

theta = OT_a.THETA/57.3;
phi = OT_a.PHI/57.3;
psi = OT_a.PSI/57.3;
save(['./data/', num2str(index),'_',num2str(DU(1)),'-',num2str(DU(2))], ...
    'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz', ...
    'fx', 'fy', 'fz', 'mx', 'my', 'mz', ...
    'u', 'v', 'w', ... 
    'u1','u2','u3','u4','v1','v2','v3','v4','w1','w2','w3','w4',...
    'omega1_bar', 'omega2_bar', 'omega3_bar', 'omega4_bar', ...
    'omega1','omega2','omega3','omega4',...
    'alpha','psi_mod','DU',...
    'p','q','r',...
    'phi','theta','psi',...
    'omega1_dot','omega2_dot','omega3_dot','omega4_dot');
end

return;

%% Project specific force f on the inertia frame
fi = zeros(size(f)); fi_model = zeros(size(f));
for i = 1:N
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
       
    R_IB = R_BI';
    fi(i,:) = (R_IB*[f(i,1:2),0]')';
    fi_model(i,:) = (R_IB*[fx(i),fy(i),0]'/m)';
end

% figure
% plot(psi_mod2(du),fi(du,1),'.'); hold on;
% plot(psi_mod2(du),fi(du,2),'.');

% figure
% subplot(2,1,1)
% plot(fi(du,1),'.'); hold on;
% plot(fi_model(du,1),'.'); 
% subplot(2,1,2)
% plot(fi(du,2),'.'); hold on;
% plot(fi_model(du,2),'.');

figure
subplot(2,1,1)
plot(-OT_a.VY_air(du),m*fi(du,1),'.'); hold on;
subplot(2,1,2)
plot(-OT_a.VY_air(du),m*fi(du,2),'.'); hold on;
%% Estimate rp by filtering the position.
Xo = butterworth(OT_a.posCO_E(:,1),4,15/526); Yo = butterworth(OT_a.posCO_E(:,2),4,15/526); Zo = butterworth(OT_a.posCO_E(:,3),4,15/526);
Xp = butterworth(Xo,4,2.5/256); Yp = butterworth(Yo,4,2.5/256); Zp = butterworth(Zo,4,2.5/256);

rp = zeros(size(VoI));
for i = DU(1):DU(2)
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
    
    R_IB = R_BI';   

    rp(i,:) = (R_BI*[Xp(i) - Xo(i), Yp(i)-Yo(i), Zp(i)-Zo(i)]')';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% below are previous codes.
%% Subtract the centrifugual force from the accelerometer measurement ( This is wrong )
rr = zeros(size(rp));
rr(:,1) = rp(:,1) - rc(1);
rr(:,2) = rp(:,2) - rc(2);
rr(:,3) = rp(:,3) - rc(3);

f = zeros(size(rp));
fI = zeros(size(rp));
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
    
    f(i,:) = (fc(i,:)' + Gr*rr(i,:)')';
    
    fI(i,:) = (R_IB*f(i,:)')';
end

F = m*f;
M = [0 0 0; diff(pqr)*512]*Iv; 

figure
plot(dVoB(du,1),'b'); hold on;plot(dVoB(du,2) ,'r');plot(dVoB(du,3) ,'g');
plot(Y(du,1),'b-.');plot(Y(du,2),'r-.');plot(Y(du,3),'g-.');
title('r_c validation');

figure
plot(zra); hold on;
plot(GGrflat*ra);
% figure
% subplot(2,1,1)
% plot(rp(du,:)); ylabel('rp')
% title('rp / rr');
% subplot(2,1,2)
% plot(rr(du,:));
% ylabel('rr');

% figure
% subplot(2,1,1)
% plot(psi_mod(du),F(du,1),'.'); hold on;
% plot(psi_mod(du),F(du,2),'.'); xlabel('\psi'); legend('Fx','Fy')
% subplot(2,1,2)
% plot(psi_mod(du),F(du,3),'.'); hold on;
% xlabel('\psi');legend('Fz');
%%
% R = zeros(size(VoI));
% for i = DU(1):DU(2)
%     theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
%     R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
%             sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
%             cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
%     
%     R_IB = R_BI';   
%     
%     R(i,:) = (R_BI*[Xo(i);Yo(i);Zo(i)])';    
% end
% dR = [0 0 0; diff(R)*512];
% ddR = [0 0 0; diff(dR)*512];
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

%% Estimate rp by analysing the accleration of point O and point P %% This is wrong!!
V_raw_inertia = [OT_a.vel_E(:,1),OT_a.vel_E(:,2),OT_a.vel_E(:,3)];
VoI = butterworth(V_raw_inertia,4,5.5/256);
VpI_filt = butterworth(VoI,4,2/256);
VoB = zeros(size(VoI));
VpB_filt = zeros(size(VpI_filt));

for i = 1:length(OT_a.TIME)
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];    
    
    R_IB = R_BI';
    
    VoB(i,:) = (R_BI * VoI(i,:)')';
    VpB_filt(i,:) = (R_BI * VpI_filt(i,:)')';
end

rpy = (VoB(:,1)-VpB_filt(:,1))./r;
rpx = (VpB_filt(:,2)-VoB(:,2))./r;
% rB = [rx, ry, zeros(size(rx))];

% rpx = mean(rpx(6000:10000));
% rpy = mean(rpy(6000:10000));
rpB = repmat([rpx, rpy, 0],length(OT_a.TIME),1);

pqr_filt = butterworth(pqr,4,15/256);
dVoI = [0,0,0;diff(VoI)]*512;
dVcI_filt = [0,0,0;diff(VpI_filt)]*512;

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
dVcB = dVoB+ cross([0,0,0;diff(pqr_filt)]*512,rB) + cross(pqr_filt,cross(pqr_filt,rpB));

% (dVc)^B/dt = d(Vc^B)/dt + omega x Vc^B; the outcome is the same as above.
VcB = VoB + cross(pqr,rB);
dVcB3 = [0,0,0;diff(VcB)]*512 + cross(pqr,VcB);

F = m*(dVcB - gB);
M = [0 0 0; diff(pqr)]*Iv;  %% warning!! the Iv is for original bb2 instead of the lightweighted one.

%% Calculating center of mass N and mm
% configure 15-6-2018
F1(1) = 9.34;
F2(1) = 193.92;
F3(1) = 187.18;
x1(1) = 0;
x2(1) = 260;
x3(1) = 0;
y1(1) = 0;
y2(1) = 0;
y3(1) = -185;

F1(2) = 193.73;
F2(2) = 4.68;
F3(2) = 193.77;
x1(2) = 0;
x2(2) = 260;
x3(2) = 260;
y1(2) = 0;
y2(2) = 0;
y3(2) = -185;

for i = 1:2
    F_total(i) = F1(i) + F2(i) + F3(i);
    x_cg(i) = (x1(i)*F1(i) + x2(i)*F2(i) + x3(i)*F3(i))/F_total(i);
    y_cg(i) = (y1(i)*F1(i) + y2(i)*F2(i) + y3(i)*F3(i))/F_total(i);    
end

%%

kp = OB_a.k_p; kq = OB_a.k_q;
p_des = OB_a.p_des; q_des = OB_a.q_des;

p_filter = butterworth(OB_a.p,4,30/256);
q_filter = butterworth(OB_a.q,4,30/256);
nup = OB_a.p_dot_des + OB_a.k_p.*(OB_a.p_des-p_filter);
nuq = OB_a.q_dot_des + OB_a.k_q.*(OB_a.q_des-q_filter);

p_dot = [0;diff(p_filter)]*512;
q_dot = [0;diff(q_filter)]*512;

figure
plot(nup); hold on;
plot(p_dot);

figure
plot(nuq); hold on;
plot(q_dot);

