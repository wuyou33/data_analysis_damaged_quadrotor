%% 
% This script estimate the forces and moments on the quadrotor by using the
% measurements from the onboard sensor. Centrifugral forces / displacement
% from c.g to the IMU are taken into account. The measuremens are stored in
% the ./data folder. The scripts in the ./system identification folder will
% use these data to identify forces and moment models.
%%
import = 0;

%%
if import 
    
    clear all;
    
    for  index = 32
        
    addpath(genpath('E:\surfdrive\DATA\_code_import_files'));
    file_index = 'E:\surfdrive\DATA\_data_OJF_june_2018\SRFDRF\Log_OJF_June';
%     file_index = 'E:\surfdrive\DATA\_data_OJF_may_2018\Log_OJF_May';
%     file_index = 'E:\surfdrive\DATA\_data_OJF_november_2017\Log_OJF_Nov';
    [OB_a,OT_a,WIND,PARA,take,DU] = import_data(index,file_index,1,'phi_ot', 1, 0, 0, 0);
    end
    
    import = 1;
end


save_force_data = 0; 
%%
DU(1)
%% select data that

% DU = [1, length(OT_a.TIME)];
% DU = [round(length(OT_a.TIME)/30), round(length(OT_a.TIME)*9.0/10)];

% DRF 2#4 removed, r>0
% DU = [25120 33700]; %#9 V = 2
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
DU = [29110 84740];

% DU = [22600 30440]; %#35 V = 2 SRF r>0 rF removed
% DU = [39420 50420]; %#35 V = 4
% DU = [52710 53420]; %#35 V = 6
% DU = [22600 53420];

% DU = [30700 42220]; %#23 v = 2 SRF r>0 lb removed May
% DU = [57670 65810]; %#23 v = 4
% DU = [77980 83230]; %#23 v = 6
% DU = [98090 100200]; %#23 v = 8
% DU = [30700 100200];

% DU = [31500 39470]; %#25 v = 2 SRF r<0 rb removed May
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
% DU = [24950 139300];

% DU = [9729 13260]; %#47 v = 1
% DU = [28250 32120]; %#47 v = 2
% DU = [43380 46690]; %#47 v = 3
% DU = [57080 60620]; %#47 v = 4
% DU = [73450 76260]; %#47 v = 5
% DU = [83860 87020]; %#47 v = 6
% DU = [93300 96690]; %#47 v = 7
% DU = [105400 106000]; %#47 v = 8
% DU = [9729 106000];

% DU = [24380 90710]; %#92 2-10 Norminal May
% DU = [24380 117810]; %#92 2-15 Norminal May

% DU = [1 63210]; %#93 4-10 Norminal beta = -90; May

% DU = [1 161100]; %#122 up_down nominal; May

% DU = [1 238300]; %#102 2-15 Norminal heavy; Nov

% DU = [8191 22960]; %#135 demo for 2018 iros paper

du = DU(1):DU(2);

%% Use EKF to estimate rc


%%
N = length(OB_a.TIME);
l = PARA.l; b = PARA.b; m = PARA.mass/1000; Iv = PARA.Iv; R = PARA.R;
Area = pi*R^2;
rho = 1.225;

omega_mean = 1200^2;

% rotor speed with respect to the body frame
rotor(1).omega_b = butterworth(OB_a.w1obs/60*2*pi,4,15/256);
rotor(2).omega_b = butterworth(OB_a.w2obs/60*2*pi,4,15/256);
rotor(3).omega_b = butterworth(OB_a.w3obs/60*2*pi,4,15/256);
rotor(4).omega_b = butterworth(OB_a.w4obs/60*2*pi,4,15/256);

rotor(1).rpm_ref = OB_a.w1ref;
rotor(2).rpm_ref = OB_a.w2ref;
rotor(3).rpm_ref = OB_a.w3ref;
rotor(4).rpm_ref = OB_a.w4ref;


% rotor speed with respect to the initial frame
for i = 1:4
    rotor(i).omega = rotor(i).omega_b - butterworth(OB_a.r,4,15/256);
    
    if mean(rotor(i).rpm_ref)<2800
        rotor(i).omega_b = 0*rotor(i).omega_b;
        rotor(i).omega = 0*rotor(i).omega_b;
    end
    
    rotor(i).omega_dot =  [0;diff(rotor(i).omega_b)]*512;
    rotor(i).omega_bar = rotor(i).omega.^2/omega_mean;
end


rates = [butterworth(OB_a.p,4,15/256),butterworth(OB_a.q,4,15/256), butterworth(OB_a.r,4,15/256)];

Va = [OT_a.U_air,OT_a.V_air,OT_a.W_air];

u = Va(:,1); 
v = Va(:,2);
w = Va(:,3);

va = sqrt(u.^2+v.^2+w.^2);

p = rates(:,1);
q = rates(:,2);
r = rates(:,3);

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

for i = 1:4
   rotor(i).SL = SL(i);
   rotor(i).SM = SM(i);
   rotor(i).SN = SN(i);
end

%% Calculate doa and dac
% dVoB = fc + gB - \alpha x rc - \omega x (\omega x rc)

pqr_dot = [0 0 0; diff(rates)*512];
dp = pqr_dot(:,1); dq = pqr_dot(:,2); dr = pqr_dot(:,3);
V_raw_inertia = [OT_a.vel_E(:,1),OT_a.vel_E(:,2),OT_a.vel_E(:,3)];
VoI = butterworth(V_raw_inertia,4,15/256);
% VoI = V_raw_inertia;
dVoI = [0 0 0; diff(VoI)*512];
g = 9.8124;

% pre-known bias of acc  based on measurement on the ground
if strcmp(take.drone,'bebop2_guido')
    dax=-0.198937;
    day=-0.118300;
    daz=-0.314115;
elseif strcmp(take.drone,'leon')
    dax = 0.3824; 
    day = 0.1153; 
    daz = -0.5213;
else
    dax = 0;
    day = -0.3;
    daz = -0.2;
end


% fc is the specific force on the IMU location
if mean(abs( butterworth(OB_a.az,4,15/256))) <= 2
    fcx = butterworth(OB_a.ax,4,15/256)*g - dax;
    fcy = butterworth(OB_a.ay,4,15/256)*g - day;
    fcz = butterworth(OB_a.az,4,15/256)*g - daz;
else
    fcx = butterworth(OB_a.ax,4,15/256) - dax;
    fcy = butterworth(OB_a.ay,4,15/256) - day;
    fcz = butterworth(OB_a.az,4,15/256) - daz;
end
fc = [fcx, fcy, fcz];

gB = zeros(size(VoI));
dVoB = zeros(size(VoI));
gI = [0 0 9.8124]';

GGr1 = zeros(size(VoI));
GGr2 = zeros(size(VoI));
GGr3 = zeros(size(VoI));
Gr = zeros(3,3,length(VoI));
R_BI_vector = zeros(3,3,N);

for i = 1:N
    theta = OT_a.THETA(i)/57.3; phi = OT_a.PHI(i)/57.3; psi = OT_a.PSI(i)/57.3;
    R_BI_vector(:,:,i) = [cos(theta)*cos(psi) cos(theta)*sin(psi) -sin(theta);
                    sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) sin(phi)*cos(theta);
                    cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi) cos(phi)*cos(theta)];                
end

for i = 1:N
    R_IB = R_BI_vector(:,:,i)';   
    R_BI = R_BI_vector(:,:,i);
    Gr1 = [0 -dr(i) dq(i); dr(i) 0 -dp(i); -dq(i) dp(i) 0];
    Gr2 = [-(q(i).^2+r(i).^2), p(i).*q(i), p(i).*r(i);
                p(i).*q(i), -(r(i).^2+p(i).^2), q(i).*r(i);
                p(i).*r(i), q(i).*r(i), -(p(i).^2+q(i).^2)];
    Gr(:,:,i) = Gr1 + Gr2;
    
    GGr1(i,:) = [0 -dr(i) dq(i)] + [-(q(i).^2+r(i).^2), p(i).*q(i), p(i).*r(i)];
    GGr2(i,:) = [dr(i) 0 -dp(i)] + [p(i).*q(i), -(r(i).^2+p(i).^2), q(i).*r(i)];
    GGr3(i,:) = [-dq(i) dp(i) 0] + [p(i).*r(i), q(i).*r(i), -(p(i).^2+q(i).^2)];
       
    gB(i,:) = (R_BI*gI)';
    
    dVoB(i,:) = R_BI*dVoI(i,:)';
 
end

GGr = [GGr1(du,:); GGr2(du,:); GGr3(du,:)];
GGrflat = [GGr1(du,:); GGr2(du,:)];
zrc = [fc(du,1)+gB(du,1)-dVoB(du,1);fc(du,2)+gB(du,2)-dVoB(du,2);fc(du,3)+gB(du,3)-dVoB(du,3)];
zra_flat = [fc(du,1); fc(du,2)];

if import == 1
   % distance from Optitrack center to imu.
    doa = lscov(GGr,zrc);

   % distance from c.g to imu
    dca = lscov(GGrflat,zra_flat);

   % validate the correctness of doa
    VcB = zeros(size(VoI));
    VoB = zeros(size(VoI));

    for i = 1:length(VoI)
        
        R_IB = R_BI_vector(:,:,i)';   
        R_BI = R_BI_vector(:,:,i);    

        Y(i,:) = (fc(i,:)' + gB(i,:)' - Gr(:,:,i)*doa)';
        VoB(i,:) = (R_BI*VoI(i,:)')';
        VcB(i,:) = VoB(i,:) + cross(rates(i,:),doa);

    end

    fig_check_doa = figure;
    plot(dVoB(du,1),'b'); hold on; plot(Y(du,1),'r-.')
    plot(dVoB(du,2) ,'b');plot(dVoB(du,3) ,'b');
    plot(Y(du,2),'r-.');plot(Y(du,3),'r-.');

end

%% plots showing the effects of doa and dca

% Compare force on c.g. and on IMU.
f = zeros(size(VoI));
for i = 1:N    
    f(i,:) = fc(i,:) + (-Gr(:,:,i)*dca)';
end

fig_force_cg = figure;
subplot(3,1,1)
plot(fc(du,1)); hold on;
plot(f(du,1)); ylabel('F_x/m [m/s^2]'); xlim([0,1500]);
 set(gca,'XTickLabel',[{'0'},{'1'},{'2'},{'3'}]);
subplot(3,1,2)
plot(fc(du,2)); hold on;
plot(f(du,2)); ylabel('F_y/m [m/s^2]'); xlim([0,1500]);
 set(gca,'XTickLabel',[{'0'},{'1'},{'2'},{'3'}]);
subplot(3,1,3)
plot(fc(du,3)); hold on;
plot(f(du,3)); ylabel('F_z [m/s^2]'); xlim([0,1500]);
 set(gca,'XTickLabel',[{'0'},{'1'},{'2'},{'3'}]);
legend('mea','corrected')
xlabel('time [s]');
set(gcf,'position',[722    27   299   306]);

% Compare velocities of rotors with / without corrections
d_real = [ l + doa(1),    -b + doa(2), 0 + doa(3);
           l + doa(1),     b + doa(2), 0 + doa(3);
          -l + doa(1),     b + doa(2), 0 + doa(3);
          -l + doa(1),    -b + doa(2), 0 + doa(3)];
 
d_fake = [ l ,    -b , 0;
           l ,     b , 0;
          -l ,     b , 0;
          -l ,    -b , 0];

V = cross(rates,repmat(doa',size(Va,1),1)) + Va;
for i = 1:4
   rotor(i).V =  cross(rates,repmat(d_real(i,:),size(Va,1),1)) + Va;
   rotor(i).u = rotor(i).V(:,1);
   rotor(i).v = rotor(i).V(:,2);
   rotor(i).w = rotor(i).V(:,3);
   rotor(i).va = sqrt(rotor(i).u.^2+rotor(i).v.^2+rotor(i).w.^2);
   
   rotor(i).V_fake =  cross(rates,repmat(d_fake(i,:),size(Va,1),1)) + Va;
   rotor(i).u_fake = rotor(i).V_fake(:,1);
   rotor(i).v_fake = rotor(i).V_fake(:,2);
   rotor(i).w_fake = rotor(i).V_fake(:,3);
end

fig_rotor_velocity = figure;
fig_rotor_velocity1 = subplot(3,1,1);
plot(rotor(1).u(du)); hold on; 
plot(rotor(1).u_fake(du));
xlim([0,1500]); ylabel('u_1 [m/s]');
set(gca,'XTickLabel',[{'0'},{'1'},{'2'},{'3'}]);
fig_rotor_velocity2 = subplot(3,1,2);
plot(rotor(1).v(du)); hold on; 
plot(rotor(1).v_fake(du));
xlim([0,1500]);  ylabel('v_1 [m/s]');
set(gca,'XTickLabel',[{'0'},{'1'},{'2'},{'3'}]);
fig_rotor_velocity3 = subplot(3,1,3); 
plot(rotor(1).w(du)); hold on; 
plot(rotor(1).w_fake(du));
xlim([0,1500]); ylabel('w_1 [m/s]');
set(gca,'XTickLabel',[{'0'},{'1'},{'2'},{'3'}]);
xlabel('time [s]');
set(gcf,'position',[722    27   299   306]);

fig_drone_velocity = figure;
fig_drone_velocity1 = subplot(3,1,1);
plot(Va(du,1)); hold on;
plot(V(du,1));
xlim([0,1500]); ylabel('u [m/s]');
set(gca,'XTickLabel',[{'0'},{'1'},{'2'},{'3'}]);
fig_drone_velocity2 = subplot(3,1,2);
plot(Va(du,2)); hold on;
plot(V(du,2));
xlim([0,1500]);  ylabel('v [m/s]');
set(gca,'XTickLabel',[{'0'},{'1'},{'2'},{'3'}]);
fig_drone_velocity3 = subplot(3,1,3);
plot(Va(du,3)); hold on;
plot(V(du,3));
xlim([0,1500]); ylabel('w [m/s]');
set(gca,'XTickLabel',[{'0'},{'1'},{'2'},{'3'}]);
xlabel('time [s]');
set(gcf,'position',[722    27   299   306]);

%% Calculate force and moment from the measurements and single rotor data

F = m*f;
M = [0 0 0; diff(rates)*512]*Iv + cross(rates',Iv*rates')' + ...
    [Ip*q.*(-rotor(1).omega_b+rotor(2).omega_b-rotor(3).omega_b+rotor(4).omega_b),...
    -Ip*p.*(-rotor(1).omega_b+rotor(2).omega_b-rotor(3).omega_b+rotor(4).omega_b),...
     zeros(size(p))];
M(:,3) = M(:,3) + Ip * (-rotor(1).omega_dot+rotor(2).omega_dot-rotor(3).omega_dot+rotor(4).omega_dot);

Fx = f(:,1)*m;
Fy = f(:,2)*m;
Fz = f(:,3)*m;
Mx = M(:,1);
My = M(:,2);
Mz = M(:,3);

load('E:\surfdrive\DATA\AeroModels\Cz_airframe_BB2');
load('E:\surfdrive\DATA\AeroModels\Ct_model_BB2_v2.mat');
load('E:\surfdrive\DATA\AeroModels\Cq_model_BB2_v2.mat');

alpha = atan(w./(sqrt(u.^2+v.^2)))*57.3;

fz_airframe = (u.^2 + v.^2 + w.^2).* interp1(AoA_airframe, Cz_airframe, alpha ,'spline');

fx = zeros(size(u));
fy = zeros(size(u));
for i = 1:4
    rotor(i).alpha = atan(rotor(i).w./sqrt(rotor(i).u.^2+rotor(i).v.^2))*57.3;
    rotor(i).vv = rotor(i).va./rotor(i).omega/R;
    rotor(i).i_inf = intersect(find(rotor(i).vv ~= inf),find(rotor(i).vv ~= -inf));
    
    rotor(i).Ct = P52CtCq(rotor(i).alpha(rotor(i).i_inf),rotor(i).vv(rotor(i).i_inf))*k_Ct0;
    rotor(i).Cq = P52CtCq(rotor(i).alpha(rotor(i).i_inf),rotor(i).vv(rotor(i).i_inf))*k_Cq0;
    
    if isempty(rotor(i).Ct)
        rotor(i).Ct = zeros(size(rotor(i).vv));
    end
    if isempty(rotor(i).Cq)
        rotor(i).Cq = zeros(size(rotor(i).vv));
    end
    rotor(i).T = rotor(i).Ct.*(rotor(i).omega*R).^2*Area*rho;    
    rotor(i).Q = rotor(i).SN * rotor(i).Cq;
end

% Ct1 = zeros(size(vv1)); Ct2 = zeros(size(vv2)); Ct3 = zeros(size(vv3)); Ct4 = zeros(size(vv4));

fz = -(rotor(1).T + rotor(2).T + rotor(3).T+ rotor(4).T);

mx = rotor(1).SL*b*rotor(1).T + rotor(2).SL*b*rotor(2).T + rotor(3).SL*b*rotor(3).T + rotor(4).SL*b*rotor(4).T;
my = rotor(1).SM*l*rotor(1).T + rotor(2).SM*l*rotor(2).T + rotor(3).SM*l*rotor(3).T + rotor(4).SM*l*rotor(4).T;
mz = rotor(1).Q + rotor(2).Q + rotor(3).Q + rotor(4).Q;

% Fxyz Mxyz are measured (reconstructed) forces/ moments from the data.
% fxyz mxyz are forces/moments from the simple single rotor model.
if save_force_data == 1
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
        'omega1_dot','omega2_dot','omega3_dot','omega4_dot',...
        'doa');
end


%% Sensitivity analysis of doa to L_bo
% 13-5-2019
fig_DeltaDoa = figure;
% plot(doa); hold on

fig_DeltaViNorm = figure;
plot([0,0.35],ones(1,2)*0.05,'k-.','linewidth',2); hold on;

fig_delta_vi = figure;
fig_delta_vi1 = subplot(3,1,1);
plot(abs(rotor(1).u(du)-rotor(1).u_fake(du))); hold on;
fig_delta_vi2 = subplot(3,1,2);
plot(abs(rotor(1).v(du)-rotor(1).v_fake(du))); hold on;
fig_delta_vi3 = subplot(3,1,3);
plot(abs(rotor(1).w(du)-rotor(1).w_fake(du))); hold on;    
for k = 1:0.25:10
    theta_bo    = k/57.3;
    phi_bo      = k/57.3; 
    psi_bo      = k/57.3;
    
%     theta_bo = theta_bo*rand(1); phi_bo = phi_bo*rand(1); psi_bo=psi_bo*rand(1);

    dL_bo =  [0       ,      phi_bo  ,      -theta_bo;
             -phi_bo ,      0       ,       phi_bo ;
             theta_bo,      -phi_bo ,       0];
    
    deltaDoaY = zeros(3,N);
    for i = 1:N
       Omega = [p(i);q(i);r(i)]; 
       L_oi =  R_BI_vector(:,:,i);
       deltaDoaY(:,i) = dL_bo * L_oi* ([0;0;g] - dVoB(i,:)');
    end
    deltaDoaZ = [deltaDoaY(1,du)';deltaDoaY(2,du)';deltaDoaY(3,du)'];
    DeltaDoa = lscov(GGr,deltaDoaZ);

    deltaVi = zeros(3,N);
    relativeErrorVi = zeros(1,N);
    for i = 1:N
       Omega = [p(i);q(i);r(i)]; 
       L_oi =  R_BI_vector(:,:,i);   
       deltaVi(:,i) = cross(Omega,DeltaDoa) + dL_bo*L_oi*VoI(i,:)';
       relativeErrorVi(i) = norm(deltaVi(:,i))/norm(rotor(1).V(i,:));
    end

    figure(fig_DeltaDoa)
%     plot(DeltaDoa,'r:')
    plot(k,norm(DeltaDoa)/norm(doa),'ko'); hold on;

    figure(fig_delta_vi)
    subplot(fig_delta_vi1)
%     plot(abs(deltaVi(1,du)),'r:');
    subplot(fig_delta_vi2)
    plot(abs(deltaVi(2,du)),'r:');
    subplot(fig_delta_vi3)
    plot(abs(deltaVi(3,du)),'r:');
    
    figure(fig_DeltaViNorm)
%     plot(norm(dL_bo),norm(deltaVi(:,du),2)/norm(rotor(1).V(du,:),2),'bo'); hold on;
    plot(norm(dL_bo),mean(relativeErrorVi),'b.','markersize',2); hold on;
end
figure(fig_DeltaViNorm)
xlabel('|\Delta L_{BI}|'); legend('max E[ |\Delta V| / |V| ]','E[ |\Delta V| / |V|] ')
ylim([0,0.1])
%% Project specific force f on the inertia frame
fi = zeros(size(f)); fi_model = zeros(size(f));
for i = 1:N
    R_IB = R_BI_vector(:,:,i)';
    R_BI = R_BI_vector(:,:,i);
    fi(i,:) = (R_IB*[f(i,1:2),0]')';
    fi_model(i,:) = (R_IB*[fx(i),fy(i),0]'/m)';
end
% 
% figure
% subplot(2,1,1)
% plot(fi(du,1),'.'); hold on;
% plot(fi_model(du,1),'.'); 
% subplot(2,1,2)
% plot(fi(du,2),'.'); hold on;
% plot(fi_model(du,2),'.');
% 
% figure
% subplot(2,1,2)
% plot(-OT_a.VY_air(du),-m*fi(du,1),'.'); hold on;
% plot(-OT_a.VY_air(du),-m*fi_model(du,1),'.'); hold on;
% subplot(2,1,1)
% plot(-OT_a.VY_air(du),m*fi(du,2),'.'); hold on;
% plot(-OT_a.VY_air(du),m*fi_model(du,2),'.'); hold on;
