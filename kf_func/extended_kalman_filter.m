function EKF = extended_kalman_filter(IMU ,OT, DU)
%% Load all data, compute useful variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filterType = 11; % 4, 10, 11
%savePlots = 0;
g       = 9.81;
r2d     = 180/pi;
d2r     = pi/180;

% General variables
dtvec   = diff(IMU.TIME);
NN       = length(IMU.TIME);
N = DU(2)-DU(1)+1;
t       = IMU.TIME;
dt      = median(dtvec);

% Velocity in body-coordinates; gyro rates in earth-coordinates (using OT attitude angles)
% NB the scripts Rot_x etc. define rotations in opposite directions to convention!
vB_OT_OTatt     = zeros(NN,3);
rates_gyro_E    = zeros(NN,3);
for kk = 1:NN
    Rot_B_to_E          = Rot_z(pi/180*OT.PSI(kk))*Rot_y(pi/180*OT.THETA(kk))*Rot_x(pi/180*OT.PHI(kk));
    Rot_E_to_B          = Rot_x(-pi/180*OT.PHI(kk))*Rot_y(-pi/180*OT.THETA(kk))*Rot_z(-pi/180*OT.PSI(kk));
    vB_OT_OTatt(kk,:)   = (Rot_E_to_B*[OT.velCG_E(kk,1); OT.velCG_E(kk,2); OT.velCG_E(kk,3)])';    
    rates_gyro_E(kk,:)  = (Rot_B_to_E*[IMU.P(kk); IMU.Q(kk); IMU.R(kk)])'; % DEG
end

% Initial attitude and velocity
att0            = pi/180*[OT.PHI(1); OT.THETA(1); OT.PSI(1)]; % initialise with optitrack (RAD)
Rot_E_to_B_0    = Rot_x(-att0(1))*Rot_y(-att0(2))*Rot_z(-att0(3)); % Initial E->B rotation matrix
v0              = Rot_E_to_B_0*[OT.velCG_E(1,1);OT.velCG_E(1,2);OT.velCG_E(1,3)]; % initial V in body coordinates

% integrated gyro measurements (DEG)
integ_IMU_omx   = integral_DA2(IMU.P,IMU.TIME,r2d*att0(1),'rk');
integ_IMU_omy   = integral_DA2(IMU.Q,IMU.TIME,r2d*att0(2),'rk'); 
integ_IMU_omz   = integral_DA2(IMU.R,IMU.TIME,r2d*att0(3),'rk'); 

% differentiate gyro measurements
% p_filt = butterworth(IMU.P,4,30/512);
% q_filt = butterworth(IMU.Q,4,30/512);
% r_filt = butterworth(IMU.R,4,30/512);
diff_IMU_omx  = derivative(IMU.P,IMU.TIME);
diff_IMU_omy  = derivative(IMU.Q,IMU.TIME);
diff_IMU_omz  = derivative(IMU.R,IMU.TIME);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% x = [phi theta psi u v w bp bq br bax bay baz]
% z = [phi theta psi u* v* w*] (optitrack)
% u = [p q r ax ay az] (gyro + accelerometer)
%

% process and measurement equations
F   = @kf_proc_full_IO;     % process equations
H   = @kf_meas_full_IO;     % measurement equations
DFX = @kf_proc_jac_full_IO; % Jacobian of process eq
DHX = @kf_meas_jac_full_IO; % Jacobian of measurement eq
GM  = @kf_g_full_IO;        % G matrix

% system properties
nx   = 15; % number of states 
ny   = 6; % number of measurements
nu   = 6; % number of inputs

% initial state
x0  = zeros(nx,1);%[att0;0;0]; % initial states
% x0 = [d2r*[OT.PHI(DU(1)) OT.THETA(DU(1)) OT.PSI(DU(1))]];
% measured input
g = 9.81;
U_k = [[IMU.P IMU.Q IMU.R] [IMU.AX IMU.AY IMU.AZ] [diff_IMU_omx diff_IMU_omy diff_IMU_omz]]; 

% initial estimate for covariance matrix
stdx_0  = 1;%10;
P_0     = stdx_0^2;

% measured output
Z_k = [d2r*[OT.PHI OT.THETA OT.PSI] vB_OT_OTatt];

% initial estimate for covariance matrix
stdx_0  = 1;%10;
P_0     = stdx_0^2;


% % noise properties

% Wa = 6;
% Wv = 20;
% Wpqr = 36;
% Watt = 50;

% U_k_filt(:,1) = butterworth(U_k(:,1),4,40/512);
% U_k_filt(:,2) = butterworth(U_k(:,2),4,40/512);
% U_k_filt(:,3) = butterworth(U_k(:,3),4,40/512);
% U_k_filt(:,4) = butterworth(U_k(:,4),4,10/512);
% U_k_filt(:,5) = butterworth(U_k(:,5),4,10/512);
% U_k_filt(:,6) = butterworth(U_k(:,6),4,10/512);
% 
% Z_k_filt(:,1) = butterworth(Z_k(:,1),4,15/512);
% Z_k_filt(:,2) = butterworth(Z_k(:,2),4,15/512);
% Z_k_filt(:,3) = butterworth(Z_k(:,3),4,15/512);
% Z_k_filt(:,4) = butterworth(Z_k(:,4),4,20/512);
% Z_k_filt(:,5) = butterworth(Z_k(:,5),4,20/512);
% Z_k_filt(:,6) = butterworth(Z_k(:,6),4,20/512);
% % 
% Q = diag(diag(cov(U_k-U_k_filt)));
% R = diag(diag(cov(Z_k-Z_k_filt)));

% Q = 5*diag(std(U_k));
% R = 2*diag(std(Z_k));
% Q = 3*diag(std(U_k(:,:)));
% % R = diag(std(Z_k(:,:)));
% %R(3,3)=R(3,3)./range(Z_k(:,3));
% R = diag(std(Z_k(1:512,:)));
% R = 0.8*R; 

Q = 10*0.2 * diag([5^2 5^2 5^2 25^2 25^2 25^2]);
R0 = diag([(d2r*0.5)^2 (d2r*0.5)^2 (d2r*0.5)^2 0.05^2 0.05^2 0.05^2]);
R = 0.01*R0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extended Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xhat_k = zeros(NN,nx);
Err_k = zeros(NN,ny);
Std_corr_k = zeros(NN,nx);
x = x0;
x_k_1k_1 = x0;  %Ex_0; % x(0|0)=E{x_0}
P_k_1k_1 = P_0;
% P_k_1k_1 = 10;
ntt = N; 

tic;
disp('Started the EKF')
for k = DU(1):DU(2)
    
    %%% Prediction 
    x_kk_1 = rk4_gen(F,0,x_k_1k_1,0,U_k(k,:),dt); % predicted state x(k+1|k)     
    z_kk_1 = H(x_kk_1,U_k(k,:),dt); % predicted output z(k+1|k)
    Fx = DFX(x_kk_1, U_k(k,:),dt);  % Jacobian of F(x,u,t)
    G = GM(x_kk_1); % G matrix
    [Phi, Gamma] = c2d(Fx, G, dt); % discretise
    P_kk_1 = Phi*P_k_1k_1*Phi' + Gamma*Q*Gamma'; % prediction cov. matrix P(k+1|k)
    
    %%% Correction
    Hx = DHX(x_kk_1,U_k(k,:),dt); % Jacobian of H(x,u,t)
    Ve = (Hx*P_kk_1 * Hx' + R); % Pz(k+1|k) (covariance matrix of innovation)
    K = P_kk_1 * Hx' / Ve; % K(k+1) (Kalman gain)
    x_k_1k_1 = x_kk_1 + K * (Z_k(k,:) - z_kk_1)'; % estimated state x(k+1|k+1) 
    P_k_1k_1 = (eye(nx) - K*Hx) * P_kk_1 * (eye(nx) - K*Hx)' + K*R*K';  
    P_corr = diag(P_k_1k_1);
    std_corr = sqrt(diag(P_k_1k_1));
    Err_k(k,:) = Z_k(k,:) - z_kk_1;
    
    % save results
    Xhat_k(k,:) = x_k_1k_1;
    Std_corr_k(k,:) = std_corr';
    
    if rem(k-DU(1),round(ntt/50)) == 0; disp(['EKF Progress: ' num2str(round((k-DU(1))/ntt*100)) '%']); end
end
disp('Finished the EKF')
time_kf = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correct sensor measurements with estimated biases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

rate_bias = Xhat_k(:,7:9);
acc_bias = Xhat_k(:,10:12);

p_c = IMU.P + r2d*rate_bias(:,1);
q_c = IMU.Q + r2d*rate_bias(:,2);
r_c = IMU.R + r2d*rate_bias(:,3);
ax_c = g*IMU.AX + acc_bias(:,1);
ay_c = g*IMU.AY + acc_bias(:,2);
az_c = g*IMU.AZ + acc_bias(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign states to standard variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EKF.PHI = OT.PHI;
EKF.THETA = OT.THETA;
EKF.PSI = OT.PSI;
EKF.U = OT.U;
EKF.V = OT.V;
EKF.W = OT.W;

EKF.PHI = Xhat_k(:,1);
EKF.THETA = Xhat_k(:,2);
EKF.PSI = Xhat_k(:,3);
EKF.U = Xhat_k(:,4);
EKF.V = Xhat_k(:,5);
EKF.W = Xhat_k(:,6);

EKF.rate_bias = zeros(1,NN);
EKF.acc_bias = zeros(1,NN);

EKF.rate_bias = rate_bias;
EKF.acc_bias = acc_bias;
EKF.ru = Xhat_k(:,13:15);
end

